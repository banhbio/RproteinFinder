const OneOrVector{T} = Union{T, Vector{T}}

function builddatabase!(; sources::OneOrVector{Tuple{String,Tuple{String,Tuple{Int,Int}}}}, hmmdir::String, ko_list::String, outdir::String, cpu::Int, fromhmmresult::Bool)
    if !isa(sources, Vector)        
        sources = [sources]
    end
    kofamoutdir = joinpath(outdir, "kofam")
    mkpath(kofamoutdir)

    kofam_results = map(sources) do source
        namae = basename(first(source))
        kofamoutdirforeach = joinpath(kofamoutdir, namae)
        mkpath(kofamoutdirforeach)
        if !fromhmmresult
            (first(source), runkofamscan!(first(source), hmmdir, ko_list, kofamoutdirforeach, cpu), last(source))
        else
            kofamout = parsehmmresult(first(source), hmmdir, ko_list, kofamoutdirforeach)
            (first(source), kofamout, last(source))
        end
    end

    build!(kofam_results, outdir)
end

function runkofamscan!(source_path::String, hmmdir::String, ko_list::String, outdir::String, cpu::Int)
    kofamscan = Kofamscan(source_path, outdir, hmmdir, ko_list, cpu)
    run(kofamscan)
    kofamout = result(kofamscan)
    return kofamout
end

function parsehmmresult(source_path::String, hmmdir::String, ko_list::String, outdir::String)
    namae = basename(source_path)
    kofamout_path = joinpath(outdir, namae * ".kofam.tblout")
    kofamout = Kofamout(kofamout_path, source_path, hmmdir, ko_list)
    return kofamout
end

function build!(kofam_results::Vector{Tuple{String,Kofamout,Tuple{String,Tuple{Int,Int}}}}, outdir::String)
    fasta_out = joinpath(outdir, "rproteins.fasta")
    taxid_table = joinpath(outdir, "rproteins.taxid")

    rm(fasta_out)
    rm(taxid_table)

    taxid_tmp = taxid_table * ".tmp"
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
        hit_ids = map(x -> id(x), kofam_hits)
        
        tmp_ids_file = joinpath(outdir, basename(source) * ".hit_id")
        o_tmp = open(tmp_ids_file, "w")
        for id in hit_ids
            write(o_tmp, "$(id)\n")
        end
        close(o_tmp)

        seqkitgrep = SeqkitGrep(source, fasta_out, tmp_ids_file)
        run(seqkitgrep)

        tmpdb_path = joinpath(outdir, basename(taxid_path) * ".db")
        db = SQLite.DB(tmpdb_path)
        BlastLCA.create!(db, taxid_path; header=false, delim="\t", accession_col=accession_col, taxid_col=taxid_col)

        for hit in hit_ids
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