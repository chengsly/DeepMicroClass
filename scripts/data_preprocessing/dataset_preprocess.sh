#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pytorch

python scripts/data_preprocessing/remove_plasmid.py -i data/raw_sequence/with_plasmid -o data/raw_sequence/without_plasmid
mv data/raw_sequence/without_plasmid/pre20_* data/clean_pre20
mv data/raw_sequence/without_plasmid/post20_* data/clean_post20

mkdir -p data/filtered_fasta

python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.1 -i1 data/clean_pre20/pre20_archaea.fa -i2 data/clean_post20/post20_archaea.fa -o data/filtered_fasta/archaea.fa -d data/mash/dist/archaea.tsv
python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.1 -i1 data/clean_pre20/bacteria_no_plasmid.fa -i2 data/clean_post20/bacteria_no_plasmid.fa -o data/filtered_fasta/bacteria.fa -d data/mash/dist/bacteria.tsv
cat data/filtered_fasta/archaea.fa data/filtered_fasta/bacteria.fa > data/filtered_fasta/prokaryote.fa

python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.1 -i1 data/clean_pre20/pre20_eukaryote.fa -i2 data/clean_post20/post20_eukaryote.fa -o data/filtered_fasta/euk.fa -d data/mash/dist/eukaryote.tsv

python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.1 -i1 data/plasmid/pre20_plasmid.fa -i2 data/plasmid/post20_plasmid.fa -o data/filtered_fasta/plasmid.fa -d data/mash/dist/plasmid.tsv
python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.1 -i1 data/clean_pre20/pre20_euk.fa -i2 data/clean_post20/post20_euk.fa -o data/filtered_fasta/euk.fa -d data/mash/dist/eukaryote.tsv

python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.1 -i1 data/vhdb/prokvir_pre20.fa -i2 data/vhdb/prokvir_post20.fa -o data/filtered_fasta/prokvir.fa -d data/mash/dist/prokvir.tsv
python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.1 -i1 data/vhdb/eukvir_pre20.fa -i2 data/vhdb/eukvir_post20.fa -o data/filtered_fasta/eukvir.fa -d data/mash/dist/eukvir.tsv

python scripts/remove_similar_seq.py --threshold 0.1 -i1 data/training_fasta/eukaryote.fa -i2 data/MMETSP_eukaryote_subsampled.fa -o data/filtered_fasta/MMETSP.fa -d data/dist/MMETSP.tsv
cat data/filtered_fasta/MMETSP.fa data/filtered_fasta/eukaryote.fa > data/filtered_fasta/euk_MMETSP.fa

dir=data/filtered_fasta_0.05

python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.05 -i1 data/clean_pre20/pre20_archaea.fa -i2 data/clean_post20/post20_archaea.fa -o ${dir}/archaea.fa -d data/mash/dist/archaea.tsv
python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.05 -i1 data/clean_pre20/bacteria_no_plasmid.fa -i2 data/clean_post20/bacteria_no_plasmid.fa -o ${dir}/bacteria.fa -d data/mash/dist/bacteria.tsv
cat ${dir}/archaea.fa ${dir}/bacteria.fa > ${dir}/prokaryote.fa

python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.05 -i1 data/clean_pre20/pre20_eukaryote.fa -i2 data/clean_post20/post20_eukaryote.fa -o ${dir}/euk.fa -d data/mash/dist/eukaryote.tsv

python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.05 -i1 data/plasmid/pre20_plasmid.fa -i2 data/plasmid/post20_plasmid.fa -o ${dir}/plasmid.fa -d data/mash/dist/plasmid.tsv
python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.05 -i1 data/clean_pre20/pre20_euk.fa -i2 data/clean_post20/post20_euk.fa -o ${dir}/euk.fa -d data/mash/dist/eukaryote.tsv

python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.05 -i1 data/vhdb/prokvir_pre20.fa -i2 data/vhdb/prokvir_post20.fa -o ${dir}/prokvir.fa -d data/mash/dist/prokvir.tsv
python scripts/data_preprocessing/remove_similar_seq.py --threshold 0.05 -i1 data/vhdb/eukvir_pre20.fa -i2 data/vhdb/eukvir_post20.fa -o ${dir}/eukvir.fa -d data/mash/dist/eukvir.tsv

python scripts/remove_similar_seq.py --threshold 0.05 -i1 data/training_fasta/eukaryote.fa -i2 data/MMETSP_eukaryote_subsampled.fa -o ${dir}/MMETSP.fa -d data/dist/MMETSP.tsv
cat ${dir}/MMETSP.fa ${dir}/eukaryote.fa > ${dir}/euk_MMETSP.fa

python scripts/data_preprocessing/split_fasta.py --input data/clean_pre20/pre20_eukaryote.fa --output_dir data/single_fasta --overwrite
python scripts/data_preprocessing/split_fasta.py --input data/filtered_fasta/euk.fa --output_dir data/single_fasta --overwrite
python scripts/data_preprocessing/split_fasta.py --input data/clean_pre20/pre20_archaea.fa --output_dir data/single_fasta --overwrite
python scripts/data_preprocessing/split_fasta.py --input data/clean_pre20/pre20_bacteria.fa --output_dir data/single_fasta --overwrite
python scripts/data_preprocessing/split_fasta.py --input data/filtered_fasta/prokaryote.fa --output_dir data/single_fasta --overwrite
python scripts/data_preprocessing/split_fasta.py --input data/plasmid/pre20_plasmid.fa --output_dir data/single_fasta --overwrite
python scripts/data_preprocessing/split_fasta.py --input data/filtered_fasta/plasmid.fa --output_dir data/single_fasta --overwrite
python scripts/data_preprocessing/split_fasta.py --input data/vhdb/prokvir_pre20.fa --output_dir data/single_fasta --overwrite
python scripts/data_preprocessing/split_fasta.py --input data/filtered_fasta/prokvir.fa --output_dir data/single_fasta --overwrite
python scripts/data_preprocessing/split_fasta.py --input data/vhdb/eukvir_pre20.fa --output_dir data/single_fasta --overwrite
python scripts/data_preprocessing/split_fasta.py --input data/filtered_fasta/eukvir.fa --output_dir data/single_fasta --overwrite


python scripts/data_preprocessing/fasta_to_onehot.py --input_dir data/single_fasta --output_dir data/single_onehot