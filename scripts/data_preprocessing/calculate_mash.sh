#!/bin/bash

working_dir=/home/tianqi/project/DeepMicrobeFinder/data
cd $working_dir

sketch_dir=/home/tianqi/project/DeepMicrobeFinder/data/mash/sketch

if [ ! -d "$sketch_dir" ]; then
    mkdir -p $sketch_dir
fi

dist_dir=/home/tianqi/project/DeepMicrobeFinder/data/mash/dist

if [ ! -d "$dist_dir" ]; then
    mkdir -p $dist_dir
fi

sketch_switch=0
dist_switch=0

if [ sketch_switch ]; then


mash sketch -o $sketch_dir/pre20_archaea -p 16 -i $working_dir/clean_pre20/pre20_archaea.fa
mash sketch -o $sketch_dir/post20_archaea -p 16 -i $working_dir/clean_post20/post20_archaea.fa

mash sketch -o $sketch_dir/pre20_bacteria -p 16 -i $working_dir/clean_pre20/pre20_bacteria.fa
mash sketch -o $sketch_dir/post20_bacteria -p 16 -i $working_dir/clean_post20/post20_bacteria.fa

mash sketch -o $sketch_dir/pre20_eukaryote -p 16 -i $working_dir/clean_pre20/pre20_euk.fa
mash sketch -o $sketch_dir/post20_eukaryote -p 16 -i $working_dir/clean_post20/post20_euk.fa

mash dist -p 16 $sketch_dir/pre20_archaea.msh $sketch_dir/post20_archaea.msh > $dist/archaea.tsv
mash dist -p 16 $sketch_dir/pre20_bacteria.msh $sketch_dir/post20_bacteria.msh > $dist/bacteria.tsv
mash dist -p 16 $sketch_dir/pre20_eukaryote.msh $sketch_dir/post20_eukaryote.msh > $dist/eukaryote.tsv