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

    @assert isfile(parsed_args["profilelist"])
    profilelist = abspath(parsed_args["profilelist"])

    @assert isfile(parsed_args["db_path"])
    db = abspath(parsed_args["db_path"])

    @assert isdir(parsed_args["taxonomy_db"])
    taxonomy_db = abspath(parsed_args["taxonomy_db"])

    @assert isfile(parsed_args["seq2taxonomy_db"])
    seq2taxonomy_db = abspath(parsed_args["seq2taxonomy_db"])

    @assert isdir(parsed_args["hmmdir"])
    hmmdir = abspath(parsed_args["hmmdir"])

    outputdir = abspath(parsed_args["outputdir"])

    minimal = 0.9
    cutoff = 0.8
    ranks =[:superkingdom, :phylum, :class, :order, :family, :genus, :species]
    precision = Dict{Symbol, Float64}(
               :class => 0.50,
               :order => 0.65,
               :family => 0.80,
               :genus => 0.95,
               :species => 1.0)

    @info "Parsed args:" thread input profilelist outputdir db taxonomy_db seq2taxonomy_db hmmdir minimal cutoff ranks presicion

    
    RproteinFinder.findrproteins(;query=input,
                                outputdir=outputdir,
                                profilelist_path=profilelist,
                                db_path=db,
                                taxonomy=taxonomy_db,
                                taxid_db=seq2taxonomy_db,
                                hmmdir=hmmdir,
                                cpu=thread,
                                blastlca_minimal=minimal,
                                blastlca_cutoff=cutoff,
                                blastlca_rank=rank,
                                blastlca_precision=precision)
end

main()
