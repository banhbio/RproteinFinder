function buildprofiles(;inputlist::String, orthoDBdir::String, outputdir::String, profilelist::String, cluster_identity::Float64, cluster_coverage::Float64, cpu::Int)

    f = open(inputlist, "r")
    mkpath(outputdir)

    for line in readlines(f)
        og_name = last(split(line, "\t"))

        fa = joinpath(orthoDBdir, og_name * ".fasta")

        records = [record for record in open(FASTA.Reader, fa)]
        filtered_records = remove_2σ(records)
        filtered = joinpath(outputdir, og_name * ".filtered.fasta")
        writer = open(FASTA.Writer, filtered)
        for record in filtered_records
            write(writer, record)
        end
        close(writer)

        reduced = joinpath(outputdir, og_name * ".reduced.fasta")

        cdhit = Cdhit(filtered, reduced, cpu, cluster_identity, cluster_coverage)
        run(cdhit)

        reduced_2 = joinpath(outputdir, og_name * ".reduced_2.fasta")
        trimal = Trimal(reduced, reduced_2)
        run(trimal)

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
    end

    close(f)
end

function makeadaptivethreshold(;inputlist::String, orthoDBdir::String, allfasta::String, outputdir::String, profilelist::String, cluster_identity::Float64, cluster_coverage::Float64, cpu::Int)
    f = open(inputlist, "r")

    allrecords = [record for record in open(FASTA.Reader, allfasta)]
    mkpath(outputdir)

    o = open(profilelist, "w")
    for line in readlines(f)
        (rprotein_name, og_name) = last(split(line, "\t"))

        fa = joinpath(orthoDBdir, og_name * ".fasta")
        true_records = [record for record in open(FASTA.Reader, fa)]
        true_records_id = map(x -> identifier(x), true_records)

        false_records = [record for record in allrecords if ! in(identifier(record), true_records_id)]
        false_records_path = joinpath(outputdir, og_name, "false_records.fasta")
        writer = open(FASTA.Writer, false_records_path)
        for record in false_records
            write(writer, record)
        end
        close(writer)

        shuffle!(true_records)

        len = length(true_records)
        section, extra = divrem(len, 3)
        splite_site = vcat(repeat([section + 1], extra), repeat([section], 3 - extra))

        first_site = splite_site[1]
        second_site = first_site + splite_site[2]

        true_records_set = [true_records[1:first_site], true_records[first_site+1:second_site], true_records[second_site+1:end]]
        
        threshold_f1score_pairs = Tuple{Float64,Float64}[]
        for i in 1:3
            round_i_tempdir = joinpath(outputdir, og_name, "round$i")

            true_records_test = true_records_set[i]
            true_records_test_path = joinpath(round_i_tempdir, "true_records_test.fasta")
            writer = open(FASTA.Writer, true_records_test_path)
            for record in true_records_test
                write(writer, record)
            end
            close(writer)

            true_records_train = vcat([true_records_set[j] for j in 1:3 if j != i]...)
            true_records_train_path = joinpath(round_i_tempdir, "true_records_train.fasta")
            writer = open(FASTA.Writer, true_records_train_path)
            for record in true_records_train
                write(writer, record)
            end
            close(writer)

            filtered_records = remove_2σ(true_records_train)
            filtered = joinpath(round_i_tempdir, "true_records_train.filtered.fasta")
            writer = open(FASTA.Writer, filtered)
            for record in filtered_records
                write(writer, record)
            end
            close(writer)

            reduced = joinpath(round_i_tempdir, "true_records_train.reduced.fasta")

            cdhit = Cdhit(filtered, reduced, cpu, cluster_identity, cluster_coverage)
            run(cdhit)

            reduced_2 = joinpath(round_i_tempdir, "true_records_train.reduced_2.fasta")
            trimal = Trimal(reduced, reduced_2)
            run(trimal)

            msa_path = joinpath(round_i_tempdir, "true_records_train.reduced.msa")
            msa = MSA(msa_path)

            muscle = Muscle(reduced, msa, cpu)
            run(muscle)

            profile_path = joinpath(round_i_tempdir, "true_records_train.reduced.hmm")
            profile = Profile(og_name, profile_path, 0)

            hmmbuild = Hmmbuild(msa, profile, cpu)
            run(hmmbuild)

            false_records_tblout_path = joinpath(round_i_tempdir,  "false_records.tbl")
            false_records_tblout = Tblout(false_records_tblout_path, false_records_path)
            hmmsearch = Hmmsearch(false_records_path, profile, false_records_tblout, cpu)
            run(hmmsearch)

            true_records_tblout_path = joinpath(round_i_tempdir,  "true_records_test.tbl")
            true_records_tblout = Tblout(true_records_tblout_path, true_records_test_path)
            hmmsearch = Hmmsearch(true_records_test_path, profile, true_records_tblout, cpu)
            run(hmmsearch)

            true_bitscores = bitscores(true_records_tblout)
            false_bitscores = bitscores(false_records_tblout)

            push!(threshold_f1score_pairs, maximize_f1score(true_bitscores, false_bitscores, length(true_records_test), length(false_records)))
        end

        mean_threshold = mean(map(x -> first(x), threshold_f1score_pairs))
        mean_f1score = mean(map(x -> last(x), threshold_f1score_pairs))
        
        write(o, "$(rprotein_name)\t$(og_name)\t$(mean_threshold)\t$(mean_f1score)\n")
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
