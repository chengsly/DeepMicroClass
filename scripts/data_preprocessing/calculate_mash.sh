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

    mash sketch -o $sketch_dir/pre20_eukaryote -p $threads -i $working_dir/clean_pre20/pre20_eukaryote.fa
    mash sketch -o $sketch_dir/post20_eukaryote -p $threads -i $working_dir/clean_post20/post20_eukaryote.fa

    mash sketch -o $sketch_dir/pre20_plasmid -p $threads -i $working_dir/plasmid/pre20_plasmid.fa
    mash sketch -o $sketch_dir/post20_plasmid -p $threads -i $working_dir/plasmid/post20_plasmid.fa

    mash sketch -o $sketch_dir/pre20_prokvir -p $threads -i $working_dir/vhdb/prokvir_pre20.fa
    mash sketch -o $sketch_dir/post20_prokvir -p $threads -i $working_dir/vhdb/prokvir_post20.fa
    mash sketch -o $sketch_dir/pre20_eukvir -p $threads -i $working_dir/vhdb/eukvir_pre20.fa
    mash sketch -o $sketch_dir/post20_eukvir -p $threads -i $working_dir/vhdb/eukvir_post20.fa

fi

mash dist -p $threads $sketch_dir/pre20_archaea.msh $sketch_dir/post20_archaea.msh -t > $dist_dir/archaea.tsv
mash dist -p $threads $sketch_dir/pre20_bacteria.msh $sketch_dir/post20_bacteria.msh -t > $dist_dir/bacteria.tsv
mash dist -p $threads $sketch_dir/pre20_eukaryote.msh $sketch_dir/post20_eukaryote.msh -t > $dist_dir/eukaryote.tsv
mash dist -p $threads $sketch_dir/pre20_plasmid.msh $sketch_dir/post20_plasmid.msh -t > $dist_dir/plasmid.tsv
mash dist -p $threads $sketch_dir/pre20_prokvir.msh $sketch_dir/post20_prokvir.msh -t > $dist_dir/prokvir.tsv
mash dist -p $threads $sketch_dir/pre20_eukvir.msh $sketch_dir/post20_eukvir.msh -t > $dist_dir/eukvir.tsv
