#!/bin/bash

input=$1
output_dir=$2
n=$3
tmp=${4:-/tmp}

mkdir -p $output_dir
if [ -d "$tmp" ]; then rm -Rf $tmp; fi
mkdir -p $tmp

dirname=$(dirname $input)
filename=$(basename -- "$input")
extension="${filename##*.}"
filename="${filename%.*}"

pyfasta split -n $n $input

find $dirname -maxdepth 1 -name "*.[0-9]*.fa" -exec mv {} $tmp \;
rm $input.flat
rm $input.gdx

run_prodigal() {
    f=$1
    fn=$(basename -- "$f")
    prodigal -i $f -a $2/${fn}.faa -o $2/${fn}.gff -p meta -f gff -q
}
export -f run_prodigal

parallel -j $n run_prodigal {} $tmp ::: $tmp/*.fa

cat $tmp/*.gff > $output_dir/contigs_genes.gff
cat $tmp/*.faa > $output_dir/contigs_proteins.faa

rm -Rf $tmp