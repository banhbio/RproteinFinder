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
        
        "--input", "-i"
            help = "input protein file (fasta format)"
            arg_type = AbstractString
            required = true

        "--taxonomy_db"
            help = "taxonomy database (directory where nodes.dmp and names.dmp)"
            arg_type = AbstractString
            required = true
        
        "--seq2taxonomy_db"
            help = "seq2taxonomy seqlite database"
            arg_type = AbstractString
            required = true

        "--profilelist"
            help = "output profile list"
            arg_type = AbstractString
            required = true

        "--hmmdir"
            help = "hmm profile directory"
            arg_type = AbstractString
            required = true

        "--outputdir"
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

    @assert isdir(parsed_args["taxonomy_db"])
    taxonomy_db = abspath(parsed_args["taxonomy_db"])

    @assert isfile(parsed_args["seq2taxonomy_db"])
    seq2taxonomy_db = abspath(parsed_args["seq2taxonomy_db"])

    @assert isdir(parsed_args["hmmdir"])
    hmmdir = abspath(parsed_args["hmmdir"])

    outputdir = abspath(parsed_args["outputdir"])

    @assert isfile(parsed_args["profilelist"])
    profilelist = abspath(parsed_args["profilelist"])

    taxonomy = Taxonomy.DB(taxonomy_db, "nodes.dmp","names.dmp")
    sqlite = SQLite.DB(seq2taxonomy_db)
    euk = Taxon(2759, taxonomy)

    @info "Parsed args:" thread input taxonomy_db seq2taxonomy_db hmmdir outputdir profilelist

    RproteinFinder.builddatabase(; source_path=input,
                                   taxonomic_scope=euk,
                                   taxonomy=taxonomy,
                                   taxid_sqlite=sqlite,
                                   profilelist_path=profilelist,
                                   hmmdir=hmmdir,
                                   outputdir=outputdir,
                                   cpu=thread)
end

main()