abstract type AbstractData end

@inline path(d::AbstractData) = d.path

struct Blastout <: AbstractData
    path::String
    query::String
    subject::String
end

@inline query(blastout::Blastout) = blastout.query
@inline subject(blastout::Blastout) = blastout.subject

struct Profile <: AbstractData
    path::String
    type::String
    threshold::String
end

struct Kofamout <: AbstractData
    path::String
    query::String
    profiledir::String
    ko_list::String
end

function hits(kofamout::Kofamout)
    open(path(kofamout), "r") do f
        result = String[]
        for l in eachline(f)
            row = split(l, "\t")
            push!(result, first(row))
        end
    end
    return result
end

struct Tblout <: AbstractData
    path::String
    query::String
    profile::Profile
end

function hits(tblout::Tblout)
    hit_list = Tuple{String, String, Float64, Float64, Float64}[]
    profile = tblout.profile
    if profile.type == "full"
        score_row = 6
    else
        score_row = 9
    end

    open(path(tblout), "r") do f
        for l in eachline(f)
            l[1] == '#' ? continue : nothing
            rows = split(l, r" +")
            id = rows[1]
            ko = rows[3]
            score = rows[score_row]
            evalue = rows[score_row-1]
            if profile.threshold >= score 
                hit = (id, ko, profile.threshold, score, evalue)
                push!(hit_list, hit)
            end
        end
    end

    return hit_list
end