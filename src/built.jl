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
        profile = Profile(profile_path, minbit=0)

        hmmbuild = Hmmbuild(msa, profile, cpu)
        run(hmmbuild)

        tblout_path = joinpath(outputdir, og_name * ".tbl")
        tblout = Tblout(tblout_path, fa)

        hmmsearch = Hmmsearch(fa, profile, tblout, cpu)
        run(hmmsearch)

        t = open(tblout, "r")

        min = 0
        for l in readlines(t)
            bitscore = split(l, "\t")[6]
            if bitscore > min
                min = bitscore
            end
        end

        min_bit = min * 0.9

        write(o, "$(rprotein_name)\t$(og_name)\t$(min_bit)\n")
    end

    close(f)
    close(o)
end