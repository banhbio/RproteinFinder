abstract type AbstractExternalProgram end

Base.run(ep::AbstractExternalProgram) = Base.run(ep.cmd)
result(ep::AbstractExternalProgram) = ep.result 


struct Kofamscan <: AbstractExternalProgram
    cmd::Cmd
    cpu::Int
    input::String
    config::String
    result::Kofamout
end

function Kofamscan(input::String, outdir::String, namae::String, profile_dir::String, ko_list::String, hmmsearch_path::String, paralell_path::String, cpu::Int)
    config = kofamconfig!(profile_dir, ko_list, hmmsearch_path, paralell_path, outdir)
    tmp_dir = joinpath(outdir, "tmp")
    mkpath(tmp_dir)
    output = joinpath(outdir, "$(namae).kofam.tblout")
    cmd = `exec_annotation -o $(output) --tmp-dir $(tmp_dir) --cpu=$(cpu) -c $(config) $(input)`
    kofamout = Kofamout(input, output, config)
    return Kofamscan(cmd, cpu, input, config, kofamout)
end

function kofamconfig!(profile_dir::String, ko_list::String, hmmsearch_path::String, paralell_path::String, outdir::String)
    config_path = joinpath(outdir, "config.yml")
    config =
    """
    profile: $(profile_dir)
    ko_list: $(ko_list)
    hmmsearch: $(hmmsearch_path)
    paralell: $(paralell_path)
    """
    open(config_path, "w") do io
        write(io, config)
    end
    return config_path
end

struct Blast <: AbstractExternalProgram
    cmd::Cmd
    cpu::Int
    input::String
    db::String
    result::Blastout
end

function Blast(input::String, db::String, output::String, evalue::Float64, cpu::Int)
    cmd = `diamond blastp --db $(db) --query $(input) --outfmt 6 --threads $(cpu) --evalue $(evalue) -k 500 --out $(output)`
    blastout = Blastout(output, input, db)
    return Blast(cmd, cpu, input, db, blastout)
end
