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
            if bitscore > min
                min = bitscore
            end
        end
        close(t)

        write(o, "$(rprotein_name)\t$(og_name)\t$(min)\n")
    end

    close(f)
    close(o)
end

function builddatabase(; source_path::String, taxonomic_scope::Taxon, id2taxid_path::String, profilelist_path::String, hmmdir::String, outputdir::String, cpu::Int)
    profilelist = profilefromlist(profilelist_path, hmmdir)

    allfasta_path = joinpath(outputdir, "rproteins.fasta")
    allfasta_writer = open(FASTA.Writer, allfasta_path)
    for profile in profilelist
        tblout_path = joinpath(outputdir, name(profile) * ".tbl")
        tblout = Tblout(tblout_path, source_path)

        hmmserach = Hmmsearch(source_path, profile, tblout, cpu)
        run(hmmserach)
        
        t = open(path(Tblout), "w")
        fasta_reader = open(FASTA.Reader, source_path)

        hit_ids = Vector{String}()
        for l in readline(t)
            hit_id = split(l, "\t")[1]
            push!(hit_ids, hit_id)
        end
        close(t)

        hits = [record for record in fasta_reader if (identifier(record) in hit_ids)]
        close(fasta_reader)
        

        euk_hits = filter(hits, taxonomic_scope)
        final_hits = remove_2σ(euk_hits)

        fasta_path = joinpath(outputdir, name(profile) * ".fasta")
        fasta_writer = open(FASTA.Writer, fasta_path)
        write(fasta_writer, final_hits)
        close(fasta_writer)

        write(allfasta_writer, final_hits)
    end
    close(allfasta_writer)
end

function profilefromlist(path::String, hmmdir::String)
    f = open(path, "r")

    l = Vector{Profile}()

    for line in readline(f)
        (rprotein_name , og_name, minbit) = split(line, "\t")
        hmm_path = joinpath(hmmdir, og_name * ".hmm")
        profile = Profile(rprotein_name, hmm_path, minbit)
        push!(l, profile)
    end

    return l
end

