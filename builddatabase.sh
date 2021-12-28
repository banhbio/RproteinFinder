#!/bin/bash
source /etc/profile.d/modules.sh

WORKDIR=/aptmp/ideas2/ban/RproteinFinder.jl

cd $WORKDIR

module load hmmer/3.3.2
module load julia/1.6.2

julia --project=@. bin/builddatabase.jl -t 64 -i --taxonomy_db --seq2taxonomy_db --hmmdir --profilelist ./data/profile_list.tsv --outputdir
