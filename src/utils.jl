function profilefromlist(path::String, hmmdir::String, threshold::Int=1)
    f = open(path, "r")

    l = Vector{Profile}()

    for line in readline(f)
        (rprotein_name , og_name, minbit) = split(line, "\t")
        hmm_path = joinpath(hmmdir, og_name * ".hmm")
        profile = Profile(rprotein_name, hmm_path, minbit*threshold)
        push!(l, profile)
    end

    return l
end

function Base.filter(v::Vector{FASTA.Record}, ts::Taxon, taxonomy::Taxonomy.DB, sqlite::SQLite.DB)
    
    filtered_v = Vector{FASTA.Record}()

    for record in v
        taxid = get(sqlite, identifier(record), nothing)
        
        if taxid === nothing
            @warn "record $(identifier(record)) has no taxid in $(sqlite.file)"
            continue
        end

        taxon = get(taxid, taxonomy, nothing)

        if taxon === nothing
            @warn "There is no taxon correspondinig to $(taxid)!"
            continue
        end

        if isdescendant(taxon, ts)
            push!(filtered_v, record)
        end
    end

    return filtered_v
end

function remove_2σ(v::Vector{FASTA.Record})
    length_v = map(x -> FASTA.seqlen(x), v)
    μ = mean(length_v)
    σ = std(length_v)

    min_length = μ - 2σ
    max_length = μ + 2σ

    filtered_v = filter(x -> (FASTA.seqlen(x) >= min_length) && (FASTA.seqlen(x) <= max_length), v)
    return filtered_v
end