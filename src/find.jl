function findrproteins(;query::String, outputdir::String, profilelist_path::String, db_path::String, taxonomy::Taxonomy.DB, taxid_db::SQLite.DB, hmmdir::String, cpu::Int, blastlca_minimal::Float64, blastlca_cutoff::Float64, blastlca_rank::Vector{Symbols}, blastlca_precision::Dict{Symbol, Float64})

    profilelist = profilefromlist(profilelist_path, hmmdir, 0.9)

    for profile in profilelist
        tblout_path = joinpath(outputdir, name(profile) * ".tbl")
        tblout = Tblout(tblout_path, query)

        hmmserach = Hmmsearch(query, profile, tblout, cpu)
        run(hmmserach)

        hits = hits(tblout)
        hits_path = joinpath(outputdir, "hits" , name(profile) * ".fasta")
        h = open(fasta.Writer, hits_path)
        write(h, hits)

        blastout_path = joinpath(outputdir, "blastout", name(profile) * ".tsv")
        blastout = Blastout(blastout_path, hits_path, db_path)

        blast = Blast(hits_path, db_path, blastout, 1e-05, cpu)
        run(blast)
        
        blastlca_path = joinpath(outputdir, "lca", name(profile) * ".tsv")
        fun = x-> weightedLCA(x, blastlca_minimal, blastlca_cutoff, blastlca_ranks, blastlca_presicion)
        
        blastLCA(;filepath=path(blastout),
                  outpath=blastlca_path,
                  sqlite=taxid_db,
                  taxonomy=taxonomy,
                  method=fun,
                  header=false,
                  ranks=blastlca_rank)
    end 
end