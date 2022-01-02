abstract type AbstractExternalProgram end

Base.run(ep::AbstractExternalProgram) = Base.run(ep.cmd)
result(ep::AbstractExternalProgram) = run(ep.result)

struct Cdhit <: AbstractExternalProgram
    cmd::Cmd
    cpu::Int
    input::String
    result::String
end

function Cdhit(input::String, result::String, cpu::Int, identity::Float64, coverage::Float64)
    cmd = `cd-hit -T $(cpu) -M 16000 -c $(identity) -aS $(coverage) -aL $(coverage) -i $(input) -o $(result)`
    return Cdhit(cmd, cpu, input, result)
end

struct Muscle <: AbstractExternalProgram
    cmd::Cmd
    cpu::Int
    input::String
    result::MSA
end

function Muscle(input::String, result::MSA, cpu::Int)
    cmd = `muscle -align $(input) -output $(path(result)) -threads $(cpu) -amino`
    return Muscle(cmd, cpu, input, result)
end

struct Hmmbuild <: AbstractExternalProgram
    cmd::Cmd
    cpu::Int
    input::MSA
    result::Profile
end

function Hmmbuild(input::MSA, result::Profile, cpu::Int)
    cmd = `hmmbuild --amino --cpu $(cpu) $(path(result)) $(path(input))`
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
    cmd = `hmmsearch --tblout $(path(result)) -T $(min) --cpu $(cpu) $(path(profile)) $(input)`
    return Hmmsearch(cmd, cpu, profile, result)
end

struct Blast <: AbstractExternalProgram
    cmd::Cmd
    cpu::Int
    input::String
    db::String
    result::Blastout
end

function Blast(input::String, db::String, output::Blastout, evalue::Float64, cpu::Int)
    cmd = `diamond blastp --db $(db) --query $(input) --outfmt 6 --threads $(cpu) --evalue $(evalue) -k0 --out $(path(output))`
    return Blast(cmd, cpu, input, db, output)
end
