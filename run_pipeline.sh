#!/usr/bin/env sh
set -e
module load /project/itmatlab/sharedmodules/use.shared
module load salmon-v1.9.0
module load STAR-v2.7.9a
module load samtools/1.11
module load htslib/1.11
module load gcc/10.2.0
module load R/4.2
source venv/bin/activate


### SETUP NOTE:
# In order to use the following --profile lsf command
# you need to follow the instructions at https://github.com/Snakemake-Profiles/lsf
# and set up LSF support for snakemake

bsub -e logs/snakemake.err \
     -o logs/snakemake.out \
     snakemake --profile lsf -j 200 -c 200 "$@"
