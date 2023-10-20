#!/bin/bash

# Process input options
while getopts "d:nt:" opt; do
    case $opt in
        d) working_dir=$OPTARG;;
        n) no_sketch=1;;
        t) threads=$OPTARG;;
        \?) echo "Invalid option -$OPTARG" >&2
            exit 1
            ;;
    esac
done

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


if [ ! $no_sketch ]; then

    mash sketch -o $sketch_dir/pre20_archaea -p $threads -i $working_dir/clean_pre20/pre20_archaea.fa
    mash sketch -o $sketch_dir/post20_archaea -p $threads -i $working_dir/clean_post20/post20_archaea.fa

    mash sketch -o $sketch_dir/pre20_bacteria -p $threads -i $working_dir/clean_pre20/pre20_bacteria.fa
    mash sketch -o $sketch_dir/post20_bacteria -p $threads -i $working_dir/clean_post20/post20_bacteria.fa

    mash sketch -o $sketch_dir/pre20_eukaryote -p $threads -i $working_dir/clean_pre20/pre20_euk.fa
    mash sketch -o $sketch_dir/post20_eukaryote -p $threads -i $working_dir/clean_post20/post20_euk.fa

    mash sketch -o $sketch_dir/euk_vir -p $threads -i $working_dir/vhdb_euk_vir.fa

fi

mash dist -p $threads $sketch_dir/pre20_archaea.msh $sketch_dir/post20_archaea.msh -t > $dist_dir/archaea.tsv
mash dist -p $threads $sketch_dir/pre20_bacteria.msh $sketch_dir/post20_bacteria.msh -t > $dist_dir/bacteria.tsv
mash dist -p $threads $sketch_dir/pre20_eukaryote.msh $sketch_dir/post20_eukaryote.msh -t > $dist_dir/eukaryote.tsv
