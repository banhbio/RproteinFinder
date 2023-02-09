using Comonicon

"""
make ribosomal protein database from inputs

# Args
- `input`: input files. should be [fasta1,taxonmap1 fasta2,taxonmap2 ...]

# Options

- `-t, --threads <int>`: is thread used in running external programs
- `--hmmdir <path>`: is the path to the hmm profile directory
- `--ko-list <path>`: is the path to the ko_list.txt file
- `--nodes-file <path>`: is the path to the nodes file
- `--names-file <path>`: is the path to the names file
- `--output <path>`: is the path to the output directory

# Flags

- `--from-hmmresult`: start from hmmsearch results
"""
@cast function build(input...;
                        hmmdir=nothing,
                        ko_list=nothing,
                        nodes_file=nothing,
                        names_file=nothing,
                        output=nothing,
                        from_hmmresult::Bool=false,
                        thread::Int=1
                    )
    @assert thread > 0

    sources = Tuple{String, String}.(split.(input, ",")) |> collect

    builddatabase!(; sources=sources,
                    hmmdir=hmmdir,
                    ko_list=ko_list,
                    outdir=output,
                    nodes=nodes_file,
                    names=names_file,
                    cpu=thread,
                    fromhmmresult=from_hmmresult
                    )
end

"""
find ribosomal protein in the input

# Args
- `input`: input fasta file.

# Options

- `-t, --threads <int>`: is thread used in running external programs
- `--hmmdir <path>`: is the path to the hmm profile directory
- `--ko-list <path>`: is the path to the ko_list.txt file
- `--nodes-file <path>`: is the path to the nodes file
- `--names-file <path>`: is the path to the names file
- `--output <path>`: is the path to the output directory
- `--tempdir <path>`: is the path to the temp directory
- `--diamond-db <path>`: is the path to the diamond database directory
"""
@cast function find(input;
                    thread::Int=1,
                    nodes_file=nothing,
                    names_file=nothing,
                    hmmdir=nothing,
                    ko_list=nothing,
                    diamond_db=nothing,
                    tempdir=nothing,
                    output=nothing
                    )
    @assert thread > 0
    input = abspath(input)
    nodes_file=abspath(nodes_file)
    names_file=abspath(names_file)
    hmmdir = abspath(hmmdir)
    ko_list = abspath(ko_list)
    diamond_db=abspath(diamond_db)
    tempdir = abspath(tempdir)
    output = abspath(output)

    RproteinFinder.findrproteins(;query=input,
                                output=output,
                                tempdir=tempdir,
                                ko_list=ko_list,
                                db_path=diamond_db,
                                hmmdir=hmmdir,
                                cpu=thread
                                )
end

"""
"""
@main
