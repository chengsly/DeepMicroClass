#!/bin/bash

# source activate def
eval "$(conda shell.bash hook)"
conda activate tf2

INPUT_DIR=$1
LENGTH=$2
ENCODE=$3

for file in $INPUT_DIR/*.fasta
do
	filename=$(basename -- "$file")
	extension="${filename##*.}"
	filename="${filename%.*}"
	echo $filename
	python encode.py -i $file -p $filename -l $LENGTH -e $ENCODE -t DNA 
	echo "finished encoding $file"
done
