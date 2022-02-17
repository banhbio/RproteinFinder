const OneOrVector{T} = Union{T, Vector{T}}

function builddatabase!(; sources::OneOrVector{Tuple{String,Tuple{String,Tuple{Int,Int}}}}, hmmdir::String, ko_list::String, outdir::String, cpu::Int)
    if !isa(sources, Vector)        
        sources = [sources]
    end

    kofamoutdir = joinpath(outdir, "kofam")
    mkpath(kofamoutdir)

    kofam_results = map(sources) do source
        (first(source), runkofamscan!(first(source), hmmdir, ko_list, kofamoutdir, cpu), last(source))
    end

    build!(kofam_results, outdir)
end

function runkofamscan!(source_path::String, hmmdir::String, ko_list::String, outdir::String, cpu::Int)
    namae = basename(source_path)
    kofamscan = Kofamscan(source_path, outdir, namae, hmmdir, ko_list, cpu)
    run(kofamscan)
    kofamout = result(kofamscan)
    return kofamout
end

function build!(kofam_results::Vector{Tuple{String,Kofamout,Tuple{String,Tuple{Int,Int}}}}, outdir::String)
    fasta_out = joinpath(outdir, "rproteins.fasta")
    taxid_table = joinpath(outdir, "rproteins.taxid")

    taxid_tmp = taxid_table * "tmp"
    o = open(taxid_tmp, "w")
    for result in kofam_results
        source = first(result)
        kofamout = result[2]
        taxid_pairs = last(result)
        taxid_path = first(taxid_pairs)
        col_pair = last(taxid_pairs)
        accession_col = first(col_pair)
        taxid_col = last(col_pair)

        kofam_hits = hits(kofamout)
        
        seqkitgrep = SeqkitGrep(source, fasta_out, kofam_hits)
        run(seqkitgrep)

        tmpdb_path = joinpath(outdir, basename(taxid_path) * ".db")
        db = SQLite.DB(tmpdb_path)
        BlastLCA.create!(db, taxid_path; header=false, delim="\t", accession_col=accession_col, taxid_col=taxid_col)

        for hit in kofam_hits
            taxid = get(db, hit, nothing)
            write(o, "$(hit)\t$(taxid)\n")
        end
        rm(tmpdb_path)
    end
    close(o)

    rm_duprow(taxid_tmp, taxid_table)

    makeblastdb = MakeBlastDB(fasta_out)
    run(makeblastdb)
    
    new_db = SQLite.DB(taxid_table * ".db")
    BlastLCA.create!(new_db, taxid_table, header=false, delim="\t", accession_col=1, taxid_col=2)
end

function rm_duprow(input::String, output::String)
    cmd = pipeline(pipeline(`cat $(input)`, `sort`, `uniq`), stdout=output)
    run(cmd)
end