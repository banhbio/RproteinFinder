abstract type AbstractData end

path(d::AbstractData) = d.path

struct Profile <: AbstractData
    name::String
    path::String
    minbit::Float64
end

name(profile::Profile) = profile.name
minbit(profile::Profile) = profile.minbit

struct MSA <: AbstractData
    path::String
end

struct Tblout <: AbstractData
    path::String
    fasta::String
end
