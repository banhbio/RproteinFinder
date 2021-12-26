abstract type AbstractExternalProgram end

run(ep::AbstractExternalProgram) = run(ep.cmd)
result(ep::AbstractExternalProgram) = run(ep.result)

struct Muscle <: AbstractExternalProgram
    cmd::Cmd
    cpu::Int
    input::String
    result::MSA
end

function Muscle(input::String, result::MSA, cpu::Int)
    cmd = `musle -align $(input) -output $(result) -threads $(cpu) -amino`
    return Musle(cmd, cpu, input, result)
end

struct Hmmbuild<: AbstractExternalProgram
    cmd::Cmd
    cpu::Int
    input::MSA
    result::Profile
end

function Hmmbuild(input::MSA, result::Profile, cpu::Int)
    cmd = `hmmbuild --amino --cpu $(cpu) $(result) $(input)`
    return Hmmbuild(cmd, cpu, input, result)
end

struct Hmmsearch <: AbstractExternalProgram
    cmd::Cmd
    cpu::Int
    input::Profile
    result::Tblout
end

function Hmmsearch(input::String, profile::Profile, result::Tblout, cpu::Int)
    min = minbit(profile)
    cmd = `hmmsearch --tblout $(result) -T $(min) --cpu $(cpu) $(profile) $(input)`
    return Hmmsearch(cmd, cpu, input, result)
end

struct Blast <: AbstractData
    cmd::Cmd
    cpu::Int
    input::String
    db::String
    result::Blastout
end

function Blast(input::String, db::String, output::Blastout, evalue::Float64, cpu::Int)
    cmd = `diamond blastp --db $(db) --query $(input) --outfmt 6 --threads $(cpu) --evalue $(evalue) --out $(output)`
    return Blast(cmd, cpu, input, db, result)
end
