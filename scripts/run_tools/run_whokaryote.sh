#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate tiara

# for f in data/sampled_sequence_subsampled/*.fa
for f in data/sampled_sequence/10000/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    whokaryote.py --contigs $f --outdir "data/result_others/result_whokaryote/${filename}" --threads 16 --minsize 1 # --prodigal_file "data/result_whokaryote/${filename}/contigs_genes.gff"
done