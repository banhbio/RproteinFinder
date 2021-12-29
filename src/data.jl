abstract type AbstractData end

@inline path(d::AbstractData) = d.path

struct Profile <: AbstractData
    name::String
    path::String
    minbit::Float64
end

@inline name(profile::Profile) = profile.name
@inline minbit(profile::Profile) = profile.minbit

struct MSA <: AbstractData
    path::String
end

struct Tblout <: AbstractData
    path::String
    fasta::String
end

@inline fasta(tblout::Tblout) = tblout.fasta

function hits(tblout::Tblout)
    t = open(path(tblout), "r")

    hit_ids = Vector{String}()
    for l in eachline(t)
        if l[1] == '#'
            continue
        end

        hit_id = split(l, r" +")[1]
        push!(hit_ids, hit_id)
    end
    close(t)

    return hit_ids
end

function minbit(tblout::Tblout)
    t = open(path(tblout), "r")

    min = Inf
    for l in eachline(t)
        if l[1] == '#'
            continue
        end

        bitscore = parse(Float64, split(l, r" +")[6])
        if bitscore < min
            min = bitscore
        end
    end
    close(t)

    return min
end
struct Blastout <: AbstractData
    path::String
    query::String
    subject::String
end

@inline query(blastout::Blastout) = blastout.query
@inline subject(blastout::Blastout) = blastout.subject
