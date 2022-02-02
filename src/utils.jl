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

function maximize_f1score(true_bitscores::Vector{Float64}, false_bitscores::Vector{Float64}, length_true::Int, length_false::Int)
    append!(true_bitscores, repeat([0], length_true - length(true_bitscores)))
    append!(false_bitscores, repeat([0], length_false - length(false_bitscores)))

    candidate_thresholds = unique(sort(vcat(true_bitscores, false_bitscores)))

    threshold_f1score_pair = Tuple{Float64, Float64}[]

    for threshold in candidate_thresholds

        tp = length(filter(x -> x > threshold, true_bitscores))
        fp = length(filter(x -> x > threshold, false_bitscores))
        fn = length(filter(x -> x < threshold, true_bitscores))

        f1score = 2*tp / (2*tp + fp + fn)

        push!(threshold_f1score_pair, (threshold,f1score))
    end

    return argmax(last,threshold_f1score_pair)
end