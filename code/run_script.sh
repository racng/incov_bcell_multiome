#!/bin/bash
# 
SNAKEDIR=/home/rng/git/snakemake-merge-wgs
CONDA_ENV=conda
mkdir $1
nohup snakemake --cores $2 \
	--directory=$SNAKEDIR \
	--snakefile=$SNAKEDIR/Snakefile \
	--use-conda --conda-prefix=$SNAKEDIR/$CONDA_ENV \
	--configfile=config_$1.yaml > $1.log & echo $! > $1.pid




