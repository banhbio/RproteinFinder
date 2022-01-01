#!/bin/bash
source /etc/profile.d/modules.sh

WORKDIR=/aptmp/ban/RproteinFinder.jl

cd $WORKDIR

module load diamond/2.0.11
module load hmmer/3.3.2
module load julia/1.6.2

julia --project=@. bin/findrproteinlca.jl -t 8 -i test/MMETSP0907.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep.rename.fasta --outputdir ./test --profilelist ./data/profile_list.tsv --db_path ./data/fasta/rproteins.fasta --taxonomy_db ./data/taxonomy --seq2taxonomy_db ./data/fasta/rproteins.taxid.db --hmmdir ./data/profile/