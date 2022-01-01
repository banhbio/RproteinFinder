function buildprofiles(;inputlist::String, orthoDBdir::String, outputdir::String, profilelist::String, cluster_identity::Float64, cluster_coverage::Float64, cpu::Int)

    f = open(inputlist, "r")
    o = open(profilelist, "w")

    mkpath(outputdir)

    for line in readlines(f)
        (rprotein_name, og_name) = split(line, "\t")

        fa = joinpath(orthoDBdir, og_name * ".fasta")

        filtered = joinpath(outputdir, og_name * ".filtered.fasta")
        records = [record for record in open(FASTA.Reader, fa)]
        filtered_records = remove_2σ(records)
        writer = open(FASTA.Writer, filtered)
        for record in filtered_records
            write(writer, record)
        end
        close(writer)

        reduced = joinpath(outputdir, og_name * ".reduced.fasta")

        cdhit = Cdhit(filtered, reduced, cpu, cluster_identity, cluster_coverage)
        run(cdhit)

        msa_path = joinpath(outputdir, og_name * ".msa")
        msa = MSA(msa_path)

        muscle = Muscle(reduced, msa, cpu)
        run(muscle)

        profile_path = joinpath(outputdir, og_name * ".hmm")
        profile = Profile(og_name, profile_path, 0)

        hmmbuild = Hmmbuild(msa, profile, cpu)
        run(hmmbuild)

        tblout_path = joinpath(outputdir, og_name * ".tbl")
        tblout = Tblout(tblout_path, reduced)

        hmmsearch = Hmmsearch(reduced, profile, tblout, cpu)
        run(hmmsearch)

        min = minbit(tblout)

        write(o, "$(rprotein_name)\t$(og_name)\t$(min)\n")
    end

    close(f)
    close(o)
end

function builddatabase(; source_path::String, taxonomic_scope::Taxon, taxonomy::Taxonomy.DB, taxid_sqlite::SQLite.DB, profilelist_path::String, hmmdir::String, outputdir::String, cpu::Int)
    profilelist = profilefromlist(profilelist_path, hmmdir, 0.9)

    tblout_dir = joinpath(outputdir, "tblout")
    mkpath(tblout_dir)

    result_dir = joinpath(outputdir, "fasta")
    mkpath(result_dir)
    result_file = joinpath(result_dir, "rproteins.tsv")
    f = open(result_file, "w")

    for profile in profilelist
        tblout_path = joinpath(tblout_dir, name(profile) * ".tbl")
        tblout = Tblout(tblout_path, source_path)

        hmmserach = Hmmsearch(source_path, profile, tblout, cpu)
        run(hmmserach)
        
        hits_id = hits(tblout)
        for hit in hits_id
            taxid = get(taxid_sqlite, hit, nothing)
        
            if taxid === nothing
                @warn "record $(hit) has no taxid in $(taxid_sqlite.file)"
                continue
            end

            taxon = get(taxid, taxonomy, nothing)

            if taxon === nothing
                @warn "There is no taxon correspondinig to $(taxid)!"
                continue
            end

            if isdescendant(taxon, taxonomic_scope)
                write(f,  "$(hit)\t$(name(profile))\n")
            end
        end
    end

    close(f)
end
