#!/bin/bash
source /etc/profile.d/modules.sh

WORKDIR=/aptmp/ideas2/ban/RproteinFinder.jl
cd $WORKDIR

module load hmmer
module load ruby

julia --project=@. bin/builddatabase.jl -t 80 --input ./source/input.json --ko_list ./source/kofam/ko_list.txt --hmmdir ./source/kofam/profiles --outdir ./db

##qsub -q cdb -l select=1:ncpus=80:mem=200gb
