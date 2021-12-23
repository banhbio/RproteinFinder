function buildprofiles(;inputlist::String, orthoDBdir::String, outputdir::String, profilelist::String, cpu::Int)

    f = open(inputlist, "r")
    o = open(profilelist, "w")


    for line in readlines(f)
        (rprotein_name, og_name) = split(line, "\t")

        fa = joinpath(orthoDBdir, og_name * ".fasta")

        msa_path = joinpath(outputdir, og_name * ".msa")
        msa = MSA(msa_path)

        musle = Musle(fa, msa, cpu)
        run(musle)

        profile_path = joinpath(outputdir, og_name * ".hmm")
        profile = Profile(og_name, profile_path, 0)

        hmmbuild = Hmmbuild(msa, profile, cpu)
        run(hmmbuild)

        tblout_path = joinpath(outputdir, og_name * ".tbl")
        tblout = Tblout(tblout_path, fa)

        hmmsearch = Hmmsearch(fa, profile, tblout, cpu)
        run(hmmsearch)

        t = open(path(tblout), "r")

        min = 0
        for l in readlines(t)
            bitscore = split(l, "\t")[6]
            if bitscore < min
                min = bitscore
            end
        end
        close(t)

        write(o, "$(rprotein_name)\t$(og_name)\t$(min)\n")
    end

    close(f)
    close(o)
end

function builddatabase(; source_path::String, taxonomic_scope::Taxon, taxonomy::Taxonomy.DB, taxid_sqlite::SQLite.DB, profilelist_path::String, hmmdir::String, outputdir::String, cpu::Int)
    profilelist = profilefromlist(profilelist_path, hmmdir, 0.9)

    allfasta_path = joinpath(outputdir, "rproteins.fasta")
    allfasta_writer = open(FASTA.Writer, allfasta_path)
    for profile in profilelist
        tblout_path = joinpath(outputdir, name(profile) * ".tbl")
        tblout = Tblout(tblout_path, source_path)

        hmmserach = Hmmsearch(source_path, profile, tblout, cpu)
        run(hmmserach)
        
        hits = hits(tblout)
        euk_hits = filter(hits, taxonomic_scope, taxonomy, taxid_sqlite)
        final_hits = remove_2σ(euk_hits)

        fasta_path = joinpath(outputdir, name(profile) * ".fasta")
        fasta_writer = open(FASTA.Writer, fasta_path)
        write(fasta_writer, final_hits)
        close(fasta_writer)

        write(allfasta_writer, final_hits)
    end
    close(allfasta_writer)    
end

function profilefromlist(path::String, hmmdir::String, threshold::Int=1)
    f = open(path, "r")

    l = Vector{Profile}()

    for line in readline(f)
        (rprotein_name , og_name, minbit) = split(line, "\t")
        hmm_path = joinpath(hmmdir, og_name * ".hmm")
        profile = Profile(rprotein_name, hmm_path, minbit*threshold)
        push!(l, profile)
    end

    return l
end

function Base.filter(v::Vector{FASTA.Record}, ts::Taxon, taxonomy::Taxonomy.DB, sqlite::SQLite.DB)
    
    filtered_v = Vector{FASTA.Record}()

    for record in v
        taxid = get(sqlite, identifier(record), nothing)
        
        if taxid === nothing
            @warn "record $(identifier(record)) has no taxid in $(sqlite.file)"
            continue
        end

        taxon = get(taxid, taxonomy, nothing)

        if taxon === nothing
            @warn "There is no taxon correspondinig to $(taxid)!"
            continue
        end

        if isdescendant(taxon, ts)
            push!(filtered_v, record)
        end
    end

    return filtered_v
end

function remove_2σ(v::Vector{FASTA.Record})
    length_v = map(v, x -> seqlen(x))
    μ = mean(length_v)
    σ = std(length_v)

    min_length = μ - 2σ
    max_length = μ + 2σ

    filtered_v = filter(v, x -> (seqlen(x) >= min_length) && (seqlen(x) <= max_length))
    return filtered_v
end