#!/bin/bash
#
#script where I will put things together, all the steps needed to quantify the metagenomic data over the tree
#
#
#

data_dir='../data/'
fgs_dir='../data/final_genomes.FGS/'
out_dir='../metaproteome_application/'


protein_extension='.fasta.FGS.faa'
bins2taxonDic_f='../data/allBin2taxon_dic.json'
allNodes2taxonDic_f='../data/allNodes2taxon_dic.json'
tree_f='../data/supermat_clean_columns_clean_rows_FastTree_internalNodesNamed_rooted.outtree'

while getopts s:o:g:t: option
do
case "${option}"
in
    s) sam_f=${OPTARG};;
    o) out_dir=${OPTARG};;
    g) genomes_dir=${OPTARG};;
    t) n_threads=${OPTARG};;
esac
done

#out_dir=$(echo $sam_f | rev | cut -d'/' -f2 | rev)'/'

echo "out_dir"
echo $out_dir
echo "out_f"
echo $out_f

python3 samRead2multiFastaGenomeAbundance.py  $sam_f $out_dir $genomes_dir

python3 get_fragment2LCA_node_assignments_report.py $sam_f $bins2taxonDic_f $allNodes2taxonDic_f $tree_f $n_threads $out_dir

python3 make_iTOL_node_mostSpecificTaxaAssignment_LCA_fragments.py $sam_f $bins2taxonDic_f $allNodes2taxonDic_f $out_dir

python3 make_iTOLSymbolFileFrom_LCA_fragments.py $sam_f $bins2taxonDic_f $allNodes2taxonDic_f $out_dir
