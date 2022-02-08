using RproteinFinder
using Taxonomy
using SQLite

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
            help = "input file .json format"
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

        "--outdir"
            help = "output directory"
            arg_type = AbstractString
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    thread = parsed_args["thread"]
    @assert thread > 0

    @assert isfile(parsed_args["input"])
    input = abspath(parsed_args["input"])

    @assert isdir(parsed_args["hmmdir"])
    hmmdir = abspath(parsed_args["hmmdir"])

    @assert isfile(parsed_args["ko_list"])
    ko_list = abspath(parsed_args["ko_list"])

    outdir = abspath(parsed_args["outdir"])

    sources = RproteinFinder.parse_input(input)

    @info "Parsed args:" input hmmdir ko_list outdir

    RproteinFinder.builddatabase!(; sources=sources,
                                     hmmdir=hmmdir,
                                     ko_list=ko_list,
                                     outdir=outdir,
                                     cpu=thread)
end

main()