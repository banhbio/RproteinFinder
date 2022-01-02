function findrproteins(;query::String, outputdir::String, profilelist_path::String, db_path::String, taxonomy::Taxonomy.DB, taxid_db::SQLite.DB, hmmdir::String, cpu::Int, blastlca_minimal::Float64, blastlca_cutoff::Float64, blastlca_ranks::Vector{Symbol}, blastlca_precision::Dict{Symbol, Float64})
    mkpath(joinpath(outputdir,"hits"))
    
    @info "Start Rproteinfinder.jl to find rproteins"
    profilelist = profilefromlist(profilelist_path, hmmdir, 0.9)

    for profile in profilelist
        @info "Starting with $profile"
        tblout_path = joinpath(outputdir, "hits", name(profile) * ".tbl")
        tblout = Tblout(tblout_path, query)

        @info "@$(profile)\tRunning hmmsearch"
        hmmserach = Hmmsearch(query, profile, tblout, cpu)
        run(hmmserach)
    end

    profile_d = Dict{String, String}()
    allhits = String[]
    for profile in profilelist
        @info "@$(profile)\tWriting hmmsearch hits to .fasta file"
        tblout_path = joinpath(outputdir, "hits", name(profile) * ".tbl")
        tblout = Tblout(tblout_path, query)
        hitid = hits(tblout)
        for hit in hitid
#            @assert !haskey(profile_d, hit)
            profile_d[hit] = name(profile)
        end
        append!(allhits, hitid)
        if isempty(hitid)
            @info "@$(profile)\tThere is no hmmsearch hit in $(fasta(tblout))"
            continue
        end
    end
    
    hits_path = joinpath(outputdir, "hits.fasta")
    writer = open(FASTA.Writer, hits_path)
    reader = open(FASTA.Reader, query)
    for record in reader
        if identifier(record) in allhits
            write(writer, record)
        end
    end
    close(writer)
    close(reader)

    @info "Running diamond blastp"
    blastout_path = joinpath(outputdir, "blastout.tsv")
    blastout = Blastout(blastout_path, hits_path, db_path)

    blast = Blast(hits_path, db_path, blastout, 1e-05, cpu)
    run(blast)
        
    blastlca_path = joinpath(outputdir, "lca.tsv")
    o = open(blastlca_path, "w")

    if filesize(path(blastout)) == 0
        @info "There is no diamond blastp hit in $(db_path)"
        return
    end

    f = open(path(blastout), "r")
    fun = x-> weightedLCA(x, blastlca_minimal, blastlca_cutoff, blastlca_ranks, blastlca_precision)
    @info "Running blastLCA"
        
    lca_ch = blastLCA(f;
              sqlite=taxid_db,
              taxonomy=taxonomy,
              method=fun,
              header=false,
              ranks=blastlca_ranks
              )

    for (qseqid, taxon, lineage) in lca_ch
        id = taxid(taxon)
        lineage_txt = sprint(io -> print_lineage(io, lineage))
        write(o, "$(qseqid)\t$(profile_d[qseqid])\t$(id)\t$(lineage_txt)\n")
    end
    close(f)
    close(o)
end