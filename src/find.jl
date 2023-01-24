function findrproteins(;query::String, output::String, tempdir::String, ko_list::String, db_path::String, taxonomy::Taxonomy.DB, taxid_db::SQLite.DB, hmmdir::String, cpu::Int, blastlca_minimal::Float64, blastlca_cutoff::Float64, blastlca_ranks::Vector{Symbol}, blastlca_precision::Dict{Symbol, Float64})
    if filesize(query) == 0
        @info "$(query) is empty"
        return
    end  
    
    @info "Start Rproteinfinder.jl to find rproteins"

    kofamoutdir = joinpath(tempdir, "kofam")   
    kofamout = runkofamscan!(query, hmmdir, ko_list, kofamoutdir, cpu)
    
    hits_path = joinpath(tempdir, "hits.fasta")
    open(FASTA.Reader, query) do reader; open(FASTA.Writer, hits_path) do o
        kofamhits = hits(kofamout)
        hit_ids = map(x -> id(x), kofamhits)
        hit_record = [record for record in reader if in(identifier(record), hit_ids)]
        map(x -> write(o, x), hit_record)
    end;end
    
    if filesize(hits_path) == 0
        @info "@all There is no hmmserach hit in $(query)"
        return
    end

    @info "Running diamond blastp"
    blastout_path = joinpath(tempdir, "blastout.tsv")
    blast = Blast(hits_path, db_path, blastout_path, 1e-05, cpu)
    run(blast)
    blastout = result(blast)

    if filesize(path(blastout)) == 0
        @info "There is no diamond blastp hit in $(db_path)"
        return
    end

    @info "running blastlca"
    fun = x-> weightedLCA(x, blastlca_minimal, blastlca_cutoff, blastlca_ranks, blastlca_precision)
       
    blastLCA(path(blastout), output;
            sqlite=taxid_db,
            taxonomy=taxonomy,
            method=fun,
            header=false,
            ranks=blastlca_ranks
    )
end

function findrproteins2(;query::String, output::String, tempdir::String, ko_list::String, db_path::String, taxonomy::Taxonomy.DB, taxid_db::SQLite.DB, hmmdir::String, cpu::Int, blastlca_ranks::Vector{Symbol})
    if filesize(query) == 0
        @info "$(query) is empty"
        return
    end  
    
    @info "Start Rproteinfinder.jl to find rproteins"

    kofamoutdir = joinpath(tempdir, "kofam")   
    kofamout = runkofamscan!(query, hmmdir, ko_list, kofamoutdir, cpu)
    
    hits_path = joinpath(tempdir, "hits.fasta")
    open(FASTA.Reader, query) do reader; open(FASTA.Writer, hits_path) do o
        kofamhits = hits(kofamout)
        hit_ids = map(x -> id(x), kofamhits)
        hit_record = [record for record in reader if in(identifier(record), hit_ids)]
        map(x -> write(o, x), hit_record)
    end;end
    
    if filesize(hits_path) == 0
        @info "@all There is no hmmserach hit in $(query)"
        return
    end

    @info "Running diamond blastp"
    blastout1_path = joinpath(tempdir, "blastout_1.tsv")
    blast1 = Blast(hits_path, db_path, blastout1_path, 1e-05, cpu)
    run(blast1)
    blastout1 = result(blast1)

    if filesize(path(blastout)) == 0
        @info "There is no diamond blastp hit in $(db_path)"
        return
    end

    @info "extract best hit"
    besthits = besthit(blastout1)
    besthit_ids = map(x -> x.qseqid, besthits)
    thresholds = Dict(map(x -> Pair(x.sseqid, x.bitscore), besthits))

    besthit_path = joinpath(tempdir, "besthit_1.fasta")
    seqkitgrep = SeqkitGrep(db_path, besthit_path, besthit_ids)
    run(seqkitgrep)

    @info "Running diamond blastp"
    blastout2_path = joinpath(tempdir, "blastout_2.tsv")
    blast2 = Blast(hits_path, db_path, blastout2_path, 1e-05, cpu)
    run(blast2)
    blastout2 = result(blast2)


    @info "running blastlca"
    open(path(blastout2), "r") do f; open(output, "w") do o
        
        fun = x -> begin
            id = first(x)
            leaves = last(x)
            threshold = thresholds[id]
            taxons = [ first(l) for l in leaves if last(l).bitscore >= threshold]
            return lca(taxons)
        end

        lca_ch = blastLCA(f;
                          sqlite=taxid_db,
                          taxonomy=taxonomy,
                          method=fun,
                          header=false,
                          ranks=blastlca_ranks,
                          rmselfhit=true
                        )

        for (qseqid, taxon, lineage) in lca_ch
            id = taxid(taxon)
            lineage_txt = sprint(io -> print_lineage(io, lineage))
            write(o, "$(qseqid)\t$(id)\t$(lineage_txt)\n")
        end
    end;end
end