const OneOrVector{T} = Union{T, Vector{T}}

function builddatabase!(; sources::OneOrVector{Tuple{String,Tuple{String,Tuple{Int,Int}}}}, hmmdir::String, ko_list::String, outdir::String, cpu::Int)
    if !isa(sources, Vector)        
        sources = [sources]
    end

    kofamoutdir = joinpath(outdir, "kofam")
    mkpath(kofamoutdir)

    kofam_results = map(sources) do source
        (first(source), runkofamscan!(first(source), hmmdir, ko_list, kofamoutdir, cpu), last(source))
    end

    build!(kofam_results, outdir)
end

function runkofamscan!(source_path::String, hmmdir::String, ko_list::String, outdir::String, cpu::Int)
    hmmsearch_path = Sys.which("hmmsearch")
    parallel_path = Sys.which("parallel")
    @assert !isnothing(hmmsearch_path)
    @assert !isnothing(parallel_path)

    namae = basename(source_path)

    kofamscan = Kofamscan(source_path, outdir, namae, hmmdir, ko_list, hmmsearch_path, parallel_path, cpu)
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

function build!(kofam_results::Vector{Tuple{String,String,Tuple{String,Tuple{Int,Int}}}}, outdir::String)
    fasta_out = joinpath(outdir, ".fasta")
    taxid_table = joinpath(outdir, ".taxid")
    open(FASTA.Writer, fasta_out) do o; open(taxid_table, "r") do p
        for result in kofam_results
            source = first(result)
            kofam_hit = result[2]
            taxid_pairs = last(result)
            taxid_path = first(taxid_pairs)
            col_pair = last(taxid_pairs)
            accession_col = first(col_pair)
            taxid_col = last(col_pair)

            open(FASTA.Reader, source) do reader ; open(kofam_hit, "r") do f; open(taxid_path , "r") do g
                hit_ids = [first(split(l, "\t")) for l in eachline(f)]
                hit_record = [record for record in reader if in(identifier(record), hit_ids)]
                for record in hit_record
                    write(o,record)
                end

                for l in eachline(g)
                    row = split(l, "\t")
                    accession = row[accession_col]
                    taxid = row[taxid_col]
                    if accession in hit_ids
                        write(p, "$(accession)\t$(taxid)\n")
                    end
                end       
            end; end; end
        end
    end; end
    db = SQLite.DB(taxid_table * ".db")
    BlastLCA.create!(db, taxid_table, header=false, delim="\t", accession_col=1, taxid_col=2)
end