#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate tf2

# cat MarchTrainingFasta/BactTr.fasta MarchTrainingFasta/BactVal.fasta > training_fasta/bacteria.fa
# cat MarchTrainingFasta/ArchaeaTr.fasta MarchTrainingFasta/archaea_val.fna > training_fasta/archaea.fa
# cat MarchTrainingFasta/ProkVirusTr.fasta MarchTrainingFasta/ProkVirusVal.fasta > training_fasta/prok_virus.fa
# cat MarchTrainingFasta/EukVirusTr.fasta MarchTrainingFasta/EukVirusVal.fasta > training_fasta/euk_virus.fa
# cat MarchTrainingFasta/BacteriaPlasmid*.fasta MarchTrainingFasta/ArchaeaPlasmid*.fasta > training_fasta/plasmid.fa
# cat MarchTrainingFasta/euk_tr_2000.fasta MarchTrainingFasta/euk_val100.fasta > training_fasta/eukaryote.fa

# mash sketch -o mash_sketch/training_bacteria -p 16 -i training_fasta/bacteria.fa
# mash sketch -o mash_sketch/training_archaea -p 16 -i training_fasta/archaea.fa
# mash sketch -o mash_sketch/training_prokvir -p 16 -i training_fasta/prok_virus.fa
# mash sketch -o mash_sketch/training_eukvir -p 16 -i training_fasta/euk_virus.fa
# mash sketch -o mash_sketch/training_plasmid -p 16 -i training_fasta/plasmid.fa
# mash sketch -o mash_sketch/training_eukaryote -p 16 -i training_fasta/eukaryote.fa

# mash sketch -o mash_sketch/testing_archaea -p 16 -i post20_archaea/post20_archaea.fa
# mash sketch -o mash_sketch/testing_bacteria -p 16 -i post20_bacteria/post20_bacteria.fa
# mash sketch -o mash_sketch/testing_eukaryote -p 16 -i post20_eukaryote/post20_euk.fa
# mash sketch -o mash_sketch/testing_plasmid -p 16 -i plasmid/post20_plasmid.fa

# mash

# mash dist -p 32 mash_sketch/training_bacteria.msh mash_sketch/testing_bacteria.msh -t > dist/bacteria.tsv
# mash dist -p 32 mash_sketch/training_archaea.msh mash_sketch/testing_archaea.msh -t > dist/archaea.tsv
# mash dist -p 32 mash_sketch/training_eukaryote.msh mash_sketch/testing_eukaryote.msh -t > dist/eukaryote.tsv
# mash dist -p 32 mash_sketch/training_plasmid.msh mash_sketch/testing_plasmid.msh -t > dist/plasmid.tsv

# mash dist -p 32 mash_sketch/training_prokvir.msh mash_sketch/vhdb_prok.msh -t > dist/prokvir.tsv
# mash dist -p 32 mash_sketch/training_eukvir.msh mash_sketch/vhdb_euk.msh -t > dist/eukvir.tsv

python scripts/remove_similar_seq.py --threshold 0.01 -i1 data/training_fasta/archaea.fa -i2 data/clean_post20/post20_archaea.fa -o data/filtered_fasta/archaea.fa -d data/dist/archaea.tsv
python scripts/remove_similar_seq.py --threshold 0.01 -i1 data/training_fasta/bacteria.fa -i2 data/clean_post20/post20_bacteria.fa -o data/filtered_fasta/bacteria.fa -d data/dist/bacteria.tsv
cat data/filtered_fasta/archaea.fa data/filtered_fasta/bacteria.fa > data/filtered_fasta/prokaryote.fa
python scripts/remove_similar_seq.py --threshold 0.01 -i1 data/training_fasta/prok_virus.fa -i2 data/org_fasta/vhdb_prok_vir.fa -o data/filtered_fasta/prok_vir.fa -d data/dist/prokvir.tsv
python scripts/remove_similar_seq.py --threshold 0.01 -i1 data/training_fasta/euk_virus.fa -i2 data/org_fasta/vhdb_euk_vir.fa -o data/filtered_fasta/euk_vir.fa -d data/dist/eukvir.tsv
python scripts/remove_similar_seq.py --threshold 0.01 -i1 data/training_fasta/plasmid.fa -i2 data/plasmid/post20_plasmid.fa -o data/filtered_fasta/plasmid.fa -d data/dist/plasmid.tsv
python scripts/remove_similar_seq.py --threshold 0.01 -i1 data/training_fasta/eukaryote.fa -i2 data/org_fasta/post20_euk.fa -o data/filtered_fasta/eukaryote.fa -d data/dist/eukaryote.tsv

# seqtk sample MMETSP_eukaryote.fa 1000000 > MMETSP_eukaryote_subsampled.fa
# mash sketch -o mash_sketch/MMETSP_subsampled -p 16 -i MMETSP_eukaryote_subsampled.fa
# mash dist -p 16 mash_sketch/training_eukaryote.msh mash_sketch/MMETSP_subsampled.msh -t > dist/MMETSP.tsv
python scripts/remove_similar_seq.py --threshold 0.01 -i1 data/training_fasta/eukaryote.fa -i2 data/MMETSP_eukaryote_subsampled.fa -o data/filtered_fasta/MMETSP.fa -d data/dist/MMETSP.tsv
cat data/filtered_fasta/MMETSP.fa data/filtered_fasta/eukaryote.fa > data/filtered_fasta/euk_MMETSP.fa

python split_fasta.py --input pre20_eukaryote/pre20_euk.fa --output_dir single_fasta
python split_fasta.py --input pre20_archaea/pre20_archaea.fa --output_dir single_fasta
python split_fasta.py --input pre20_bacteria/pre20_bacteria.fa --output_dir single_fasta

python fasta_to_onehot.py --input_dir single_fasta --output_dir single_onehot