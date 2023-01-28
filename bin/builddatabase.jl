using RproteinFinder
using Taxonomy

using Logging
using ArgParse: ArgParseSettings, parse_args, @add_arg_table

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--thread", "-t"
            help = "used thread in running external programs"
            arg_type = Int
            default = 1
        
        "--input"
            help = "input file (--input fasta1,taxonmap1 fasta2,taxonmap2 ...)"
            arg_type = AbstractString
            required = true

        "--hmmdir"
            help = "hmm profile directory"
            arg_type = AbstractString
            required = true

        "--ko_list"
            help = "ko_list.txt file"
            arg_type = AbstractString
            required = true

        "--taxonomy_db"
            help = "taxonomy database (directory where nodes.dmp and names.dmp)"
            arg_type = AbstractString
            required = true

        "--outdir"
            help = "output directory"
            arg_type = AbstractString
            required = true

        "--fromhmmresult"
            help = "start from hmmsearch result"
            action = :store_true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    thread = parsed_args["thread"]
    @assert thread > 0

    input = parsed_args["input"]

    @assert isdir(parsed_args["hmmdir"])
    hmmdir = abspath(parsed_args["hmmdir"])

    @assert isfile(parsed_args["ko_list"])
    ko_list = abspath(parsed_args["ko_list"])

    @assert isdir(realpath(parsed_args["taxonomy_db"]))
    taxonomy_db = abspath(parsed_args["taxonomy_db"])

    nodes, names = joinpath.(taxonomy_db, ["nodes.dmp","names.dmp"])

    outdir = abspath(parsed_args["outdir"])

    fromhmmresult = parsed_args["fromhmmresult"]

    sources = Tuple{String, String}.(split(split(input), ","))

    @info "Parsed args:" input hmmdir ko_list outdir

    RproteinFinder.builddatabase!(; sources=sources,
                                     hmmdir=hmmdir,
                                     ko_list=ko_list,
                                     outdir=outdir,
                                     nodes=nodes,
                                     names=names,
                                     cpu=thread,
                                     fromhmmresult=fromhmmresult)
end

main()