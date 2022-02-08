function parse_input(s::String)
    f = open(s, "r")
    d = JSON.parse(f)
    result = Tuple{String,Tuple{String,Tuple{Int,Int}}}[]
    for key in keys(d)
        two = first(d[key])
        three = last(two)
        col_pair = (three["accession"], three["taxid"])
        taxid_pair = (first(two), col_pair)
        input_pair = (key, taxid_pair)
        push!(result, input_pair)
    end
    return result
end
