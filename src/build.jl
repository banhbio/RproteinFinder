const OneOrVector{T} = Union{T, Vector{T}}

function builddatabase!(; sources::OneOrVector{Tuple{String, String}}, hmmdir::String, ko_list::String, outdir::String, nodes::String, names::String, cpu::Int, fromhmmresult::Bool)
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

    build!(kofam_results, outdir, nodes, names, cpu)
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

function build!(kofam_results::Vector{Tuple{String, Kofamout, String}}, outdir::String, nodes::String, names::String, cpu::Int)
    fasta_out = joinpath(outdir, "rproteins.fasta")
    taxid_table = joinpath(outdir, "rproteins.taxid")

    rm(fasta_out, force=true)
    rm(taxid_table, force=true)

    open(taxid_table, "w") do o
        write(o, "accession.version\ttaxid\n")
    end

    for result in kofam_results
        source = first(result)
        kofamout = result[2]
        taxid_path = last(result)

        kofam_hits = hits(kofamout)
        hit_ids = id.(kofam_hits)
       
        tmp_ids_file = joinpath(outdir, basename(source) * ".hit_id")
        open(tmp_ids_file, "w") do t
            for id in hit_ids
                write(t, "$(id)\n")
            end
        end

        seqkitgrep = SeqkitGrep(source, fasta_out, tmp_ids_file, cpu)
        run(seqkitgrep)

        cmd =pipeline(`cat $(taxid_path)`, stdout=taxid_table, append=true)
        run(cmd)
    end

    makeblastdb = MakeBlastDB(fasta_out, taxid_table, nodes, names)
    run(makeblastdb)
end
