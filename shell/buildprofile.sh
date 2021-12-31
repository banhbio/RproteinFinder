#!/bin/bash
source /etc/profile.d/modules.sh

WORKDIR=/aptmp/ideas2/ban/RproteinFinder.jl

cd $WORKDIR

module load cd-hit/4.8.1.long
module load muscle/5.0.1428
module load hmmer/3.3.2
module load julia/1.6.2

julia --project=@. bin/buildprofile.jl -t 64 -i ./data/og.tsv --orthoDBdir ./data/fasta --outputdir ./data/profile --profilelist ./data/profile_list.tsv --identity 0.8 --coverage 0.9
