#!/bin/bash


# --------- #
# FUNCTIONS #
# --------- #

# help
usage()
{
    echo -e "\nusage: $0 -i <input_assembly> -o <output_dir>\n"
}

# Kaiju
kaiju_classify(){
    # arguments
    input_assembly=$1
    KAIJU_RESULTS=$2
    prefix=$3
    threads=$4    

    # database
    nr_euk=/home/db/kaiju_db/nr_euk/kaiju_db_nr_euk.fmi
    nodes_dmp=/home/db/kaiju_db/nodes.dmp
    names_dmp=/home/db/kaiju_db/names.dmp

    # search
    kaiju -z ${threads} \
        -t ${nodes_dmp} \
        -f ${nr_euk} \
        -e 5 \
        -E 0.05 \
        -x \
        -i ${input_assembly} \
        -o ${KAIJU_RESULTS}/${prefix}_kaiju.out

    # to krona
    kaiju2krona -t ${nodes_dmp} -n ${names_dmp} -i ${KAIJU_RESULTS}/${prefix}_kaiju.out -o ${KAIJU_RESULTS}/${prefix}_kaiju.out.krona
    ktImportText -o ${KAIJU_RESULTS}/${prefix}_kaiju.out.html ${KAIJU_RESULTS}/${prefix}_kaiju.out.krona

    # summarize classification
    for rank in phylum class order family genus species; do
        kaiju2table -p -t ${nodes_dmp} -n ${names_dmp} -r ${rank} -o ${KAIJU_RESULTS}/${prefix}_kaiju_summary_${rank}.tsv ${KAIJU_RESULTS}/${prefix}_kaiju.out
    done

    # addTaxonNames
    kaiju-addTaxonNames -p -t ${nodes_dmp} -n ${names_dmp} -i ${KAIJU_RESULTS}/${prefix}_kaiju.out -o ${KAIJU_RESULTS}/${prefix}_kaiju.names.out
}


# MetaEuk
metaeuk_classify(){
    # arguments
    input_assembly=$1
    MetaEuk_RESULTS=$2
    prefix=$3
    threads=$4

    # database
    REF=/home/db/metaeuk_db/MMETSP_zenodo_3247846_uniclust90_2018_08_seed_valid_taxids
    metaeuk=/home/software/src/metaeuk/build/bin/metaeuk
    mmseqs=/home/software/src/MMseqs2/build/bin/mmseqs

    # create query db
    $mmseqs createdb ${input_assembly} ${MetaEuk_RESULTS}/queryContigDB

    # predict cds 
    $metaeuk easy-predict ${MetaEuk_RESULTS}/queryContigDB ${REF} \
        ${MetaEuk_RESULTS}/${prefix}.fa tempFolder \
        -s 4 --threads 20

    # predict taxa
    $metaeuk taxtocontig ${MetaEuk_RESULTS}/queryContigDB \
        ${MetaEuk_RESULTS}/${prefix}.fa ${MetaEuk_RESULTS}/${prefix}.fa.headersMap.tsv \
        ${REF} ${MetaEuk_RESULTS}/metaeuk_tax_results tempFolder \
        --majority 0.5 --tax-lineage --lca-mode 2 --threads $threads

    # remove query db
    mv ${MetaEuk_RESULTS}/queryContigDB.lookup ${MetaEuk_RESULTS}/backup_queryContigDB.lookup
    rm ${MetaEuk_RESULTS}/queryContigDB*
}


# ---- #
# MAIN #
# ---- #

if [ "$#" != "4" ]; then
    echo -e "\nERROR: Not enough arguments provided!"
    usage
    exit 1
fi

while [[ $# -gt 0 ]]; do
    opt="$1"
    case "$opt" in
        -i|--input_assembly )
            shift
            input_assembly=$1
            ;;
        -o | --output_dir )
            shift
            output_dir=$1
            ;;
        -h | --help )
            usage
            exit
            ;;
        * )
            usage
            exit 1
            ;;
    esac
    shift
done

# output handeling
echo "The input assembly is: $input_assembly"
echo "The output directory is: $output_dir"
mkdir -p $output_dir

# run MetaEuk
MetaEuk_RESULTS=${output_dir}/MetaEuk
mkdir -p ${MetaEuk_RESULTS}
metaeuk_classify $input_assembly $MetaEuk_RESULTS "MetaEuk_" 20

# run kaiju # 
KAIJU_RESULTS=${output_dir}/kaiju
mkdir -p ${KAIJU_RESULTS}
kaiju_classify $input_assembly $KAIJU_RESULTS "Kaiju_" 20

