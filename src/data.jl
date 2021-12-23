abstract type AbstractData end

@inline path(d::AbstractData) = d.path

struct Profile <: AbstractData
    name::String
    path::String
    minbit::Float64
end

name(profile::Profile) = profile.name
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
    t = open(path(tblout), "w")
    fasta_reader = open(FASTA.Reader, fasta(Tblout))

    hit_ids = Vector{String}()
    for l in readline(t)
        hit_id = split(l, "\t")[1]
        push!(hit_ids, hit_id)
    end

    hits = [record for record in fasta_reader if (identifier(record) in hit_ids)]
    close(t)
    close(fasta_reader)

    return hits
end