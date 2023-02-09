function findrproteins(;query::String, output::String, tempdir::String, ko_list::String, db_path::String, hmmdir::String, cpu::Int)
    if filesize(query) == 0
        @info "$(query) is empty"
        return
    end  
    
    @info "Start Rproteinfinder.jl to find rproteins"

    kofamoutdir = joinpath(tempdir, "kofam")   
    kofamout = runkofamscan!(query, hmmdir, ko_list, kofamoutdir, cpu)
    
    hits_path = joinpath(tempdir, "hits.fasta")
    FASTAReader(open(query, "r")) do reader; FASTAWriter(open(hits_path, "w")) do o
        kofamhits = hits(kofamout)
        hit_ids = id.(kofamhits)
        hit_record = [record for record in reader if in(identifier(record), hit_ids)]
        map(x -> write(o, x), hit_record)
    end;end
    
    if filesize(hits_path) == 0
        @info "@all There is no hmmserach hit in $(query)"
        return
    end

    @info "Running diamond blastp"
    blast = Blast(hits_path, db_path, output, 1e-05, cpu)
    run(blast)
    blastout = result(blast)
end
