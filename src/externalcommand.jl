abstract type AbstractExternalProgram end

Base.run(ep::AbstractExternalProgram) = Base.run(ep.cmd)
result(ep::AbstractExternalProgram) = ep.result 

struct Hmmsearch <: AbstractExternalProgram
    cmd::Base.AbstractCmd
    cpu::Int
    input::String
    profile::Profile
    result::Tblout
end

function Hmmsearch(input::String, profile::Profile, output::String, evalue::Float64, cpu::Int)
    cmd = `hmmsearch --cpu $(cpu) --tblout $(output) -E $(evalue) $(path(profile)) $(input)`
    tblout = Tblout(output, input, profile)
    return Hmmsearch(cmd, cpu, input, profile, tblout)
end

struct SeqkitGrep <: AbstractExternalProgram
    cmd::Base.AbstractCmd
    input::String
    ids::String
    result::String
    cpu::Int
end

function SeqkitGrep(input::String, output::String, ids::String, cpu::Int; append=true)
    cmd = pipeline(`seqkit grep -j $(cpu) -f $(ids) $(input)`, stdout=output, append = append)
    return SeqkitGrep(cmd, input, ids, output, cpu)
end

struct MakeBlastDB <: AbstractExternalProgram
    cmd::Base.AbstractCmd
    input::String
    result::String
end

function MakeBlastDB(input::String)
    cmd = `diamond makedb --in $(input) -d $(input)`
    result = input * ".db"
    return MakeBlastDB(cmd, input, result)
end

struct Blast <: AbstractExternalProgram
    cmd::Base.AbstractCmd
    cpu::Int
    input::String
    db::String
    result::Blastout
end

function Blast(input::String, db::String, output::String, evalue::Float64, cpu::Int)
    cmd = `diamond blastp --db $(db) --query $(input) --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids --threads $(cpu) --evalue $(evalue) -k 500 --out $(output)`
    blastout = Blastout(output, input, db)
    return Blast(cmd, cpu, input, db, blastout)
end