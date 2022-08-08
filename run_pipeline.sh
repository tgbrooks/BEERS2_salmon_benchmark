#!/usr/bin/env sh
set -e
module load /project/itmatlab/sharedmodules/use.shared
module load salmon-v1.4.0
module load STAR-v2.7.9a
module load samtools/1.11
module load htslib/1.11
source venv/bin/activate


### SETUP NOTE:
# In order to use the following --profile lsf command
# you need to follow the instructions at https://github.com/Snakemake-Profiles/lsf
# and set up LSF support for snakemake

bsub -e logs/snakemake.err \
     -o logs/snakemake.out \
     snakemake --profile lsf -j 50 -c 50 --resources beers=1 "$@"
