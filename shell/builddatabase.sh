#!/bin/bash
source /etc/profile.d/modules.sh

WORKDIR=/aptmp/ban/RproteinFinder.jl

cd $WORKDIR

module load hmmer/3.3.2
module load julia/1.6.2
module load seqkit/2.1.0

julia --project=@. bin/builddatabase.jl -t 8 -i ./data/uniref/uniref90 --taxonomy_db ./data/taxonomy --seq2taxonomy_db ./data/uniref/uniref90.taxid.db --hmmdir ./data/profile --profilelist ./data/profile_list.tsv --outputdir ./data/
cat ./data/fasta/rproteins.tsv | cut -f1 |  seqkit grep -f - ./data/uniref/uniref90 > ./data/fasta/rproteins.fasta
