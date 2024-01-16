#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate vs2

for f in data/sampled_sequence/10000/*.fa
# for f in /home/shengwei/Viromes/04_read_mapping_per_sample/assemblies_per_sample/*.fasta
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    virsorter run -w data/result_others/result_virsorter2/${filename} -i $f -j 16 --provirus-off --max-orf-per-seq 10 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae all
    # virsorter run -w data/virome/virome_virsorter2/${filename} -i $f -j 16 --provirus-off --max-orf-per-seq 10 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae all
done