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
    name::String
    type::String
    threshold::Float64
end

struct Kofamout <: AbstractData
    path::String
    query::String
    hmmdir::String
    ko_list::String
end

struct KofamResult
    id::String
    ko::String
    score::Float64
    threshold::Float64
    evalue::Float64
end

id(kofamresult::KofamResult) = kofamresult.id
Base.join(kofamresult::KofamResult, delim::AbstractString) = join([kofamresult.id, kofamresult.ko, kofamresult.score, kofamresult.threshold, kofamresult.evalue], delim)

function hits(kofamout::Kofamout)
    result = KofamResult[]
    open(path(kofamout), "r") do f
        for l in eachline(f)
            row = split(l, "\t")
            id = row[1]
            ko = row[2]
            score = parse(Float64, row[3])
            threshold = parse(Float64, row[4])
            evalue = parse(Float64, row[5])
            kofamresult = KofamResult(id, ko, score, threshold, evalue)
            push!(result, kofamresult)
        end
    end
    return result
end

struct Tblout <: AbstractData
    path::String
    query::String
    profile_list::Vector{Profile}
end

function kofamhits(tblout::Tblout)
    hit_list = KofamResult[]

    profile_dict = map(tblout.profile_list) do profile
        Pair(profile.name, profile)
    end |> Dict

    open(path(tblout), "r") do f 

        for l in eachline(f)
            l[1] == '#' ? continue : nothing
            rows = split(l, r" +")

            id = rows[1]
            ko = rows[3]
        
            profile = profile_dict[ko]

            if profile.type == "full"
                score_row = 6
            else
                score_row = 9
            end
        
            score = parse(Float64, rows[score_row])
            evalue = parse(Float64, rows[score_row-1])
            if score >= profile.threshold
                hit = KofamResult(id, ko, score, profile.threshold, evalue)
                push!(hit_list, hit)
            end
        end
    end

    return hit_list
end
