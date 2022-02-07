abstract type AbstractData end

@inline path(d::AbstractData) = d.path

struct Blastout <: AbstractData
    path::String
    query::String
    subject::String
end

@inline query(blastout::Blastout) = blastout.query
@inline subject(blastout::Blastout) = blastout.subject

struct Kofamout <:AbstractData
    path::String
    query::String
    config::String
end

@inline query(kofamout::Kofamout) = kofamout.query

function hits(kofamout::Kofamout)
    hit_list = Tuple{String,String}[]
    f = open(path(kofamout), "r")
    for l in eachline(f)
        l[1] != '*' ? continue : nothing
        row = split(l, r" +")
        push!(hit_list, (row[2], row[3]))
    end
    return hit_list
end