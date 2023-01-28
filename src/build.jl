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

    taxid_tmp = taxid_table * ".tmp"
    for result in kofam_results
        source = first(result)
        kofamout = result[2]
        taxid_path = last(result)

        kofam_hits = hits(kofamout)
        hit_ids = map(x -> id(x), kofam_hits)
        
        tmp_ids_file = joinpath(outdir, basename(source) * ".hit_id")
        open(tmp_ids_file, "w") do o
            for id in hit_ids
                write(o, "$(id)\n")
            end
        end

        seqkitgrep = SeqkitGrep(source, fasta_out, tmp_ids_file, cpu)
        run(seqkitgrep)

        open(taxid_path, "r") do f; open(taxid_tmp, "w") do o
            for line in eachline(f)
                id, taxid = split(line, "\t") .|> String
                if id in hit_ids
                    write(o, "$(id)\t$(taxid)\n")
                end
            end
        end; end
    end

    rm_duprow(taxid_tmp, taxid_table)

    makeblastdb = MakeBlastDB(fasta_out, taxid_table, nodes, names)
    run(makeblastdb)
end

function rm_duprow(input::String, output::String)
    run(pipeline(`echo "accession.version\ttaxid\n"`, stdout=output))
    cmd = pipeline(pipeline(`cat $(input)`, `sort`, `uniq`), stdout=output, append=true)
    run(cmd)
end