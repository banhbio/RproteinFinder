function profilefromlist(path::String, hmmdir::String, threshold::Float64=1)
    f = open(path, "r")

    l = Vector{Profile}()

    for line in readlines(f)
        (rprotein_name , og_name, minbit) = split(line, "\t")
        hmm_path = joinpath(hmmdir, og_name * ".hmm")
        profile = Profile(rprotein_name, hmm_path, parse(Float64, minbit)*threshold)
        push!(l, profile)
    end

    return l
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