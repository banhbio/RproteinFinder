const OneOrVector{T} = Union{T, Vector{T}}

function builddatabase!(; sources::OneOrVector{Tuple{String,Tuple{String,Tuple{Int,Int}}}}, hmmdir::String, ko_list::String, outputdir::String, cpu::Int)
    if !isa(source_path, Vector)
        sources = [sources]
    end

    kofamoutdir = joinpath(outputdir, "kofam")
    mkpath(kofamoutdir)    
    resultdir = joinpath(outputdir, "fasta")
    mkpath(resultdir)

    kofam_results = map(sources) do source
        (source, runkofamscan!(first(source), hmmdir, ko_list, kofamoutdir, cpu), last(source))
    end



end

function runkofamscan!(source_path::String, hmmdir::String, ko_list::String, outdir::String, cpu::Int)
    hmmsearch_path = Sys.which("hmmsearch")
    paralell_path = Sys.which("paralell")

    namae = basename(input)

    kofamscan = Kofamscan(source_path, outdir, namae, hmmdir, ko_list, hmmsearch_path, paralell_path, cpu)
    run(kofamscan)
    kofamout = result(kofamscan)

    hit_list = hits(kofamout)
    kofam_hit = joinpath(kofamoutdir, "$(namae).ko.txt")
    open(kofam_hit ,"r") do o
        for hit in hit_list
            write(o, "$(first(hit))\t$(last(hit))\n")
        end
    end
    return kofam_hit
end

function (kofam_results::Vector{Tuple{String,String,Tuple{String,Tuple{Int,Int}}}}, outdir::String)
    for result in kofam_results
        source = first(result)
        kofam_hit = result[2]
        taxid_pairs = last(result)

        open
end