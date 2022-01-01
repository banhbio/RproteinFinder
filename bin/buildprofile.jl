using RproteinFinder

using Logging
using ArgParse: ArgParseSettings, parse_args, @add_arg_table

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--thread", "-t"
            help = "used thread in running external programs"
            arg_type = Int
            default = 1
        
        "--input", "-i"
            help = "input og list file (tsv format)"
            arg_type = AbstractString
            required = true
        
        "--orthoDBdir"
            help = "orthoDB directory"
            arg_type = AbstractString
            required = true

        "--outputdir"
            help = "output directory"
            arg_type = AbstractString
            required = true
        
        "--profilelist"
            help = "output profile list"
            arg_type = AbstractString
            required = true

        "--identity"
            help = "identity used in cd-hit"
            arg_type = Float64
            required = true

        "--coverage"
            help = "coverage used in cd-hit"
            arg_type = Float64
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

    @assert isdir(parsed_args["orthoDBdir"])
    orthoDBdir = abspath(parsed_args["orthoDBdir"])

    outputdir = abspath(parsed_args["outputdir"])
    profilelist = abspath(parsed_args["profilelist"])

    identity = parsed_args["identity"]
    @assert 0 < identity && identity < 1

    coverage = parsed_args["coverage"]
    @assert 0 < coverage && coverage < 1

    @info "Parsed args:" thread input orthoDBdir outputdir profilelist identity coverage

    RproteinFinder.buildprofiles(; inputlist=input,
                                   orthoDBdir=orthoDBdir,
                                   outputdir=outputdir,
                                   profilelist=profilelist,
                                   cpu=thread,
                                   cluster_identity=identity,
                                   cluster_coverage=coverage)
end

main()
