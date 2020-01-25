#!/bin/bash
#
#script where I will put things together, all the steps needed to quantify the metaproteomic data over the tree
#

data_dir='../data/'
fgs_dir='../data/final_genomes.FGS/'
out_dir='../metaproteome_application/'

protein_extension='.fasta.FGS.faa'
bins2taxonDic_f='../data/allBin2taxon_dic.json'
allNodes2taxonDic_f='../data/allNodes2taxon_dic.json'
tree_f='../data/supermat_clean_columns_clean_rows_FastTree_internalNodesNamed_rooted.outtree'

while getopts p:o:t:l: option
do
case "${option}"
in
    p) peptide_seq_f=${OPTARG};;
    o) out_dir=${OPTARG};;
    t) n_cores=${OPTARG};;
    l) taxonomy_level=${OPTARG};;
esac
done

python3 proteins2genomes.py $fgs_dir $protein_extension $out_dir

python3 peptide2sequenceStringMatching_parallel.py $peptide_seq_f $fgs_dir $out_dir'protein2peptideMatches/' $n_cores

python3 peptide2genomeMapping2LCA_stringMatches.py $out_dir'peptide2allProteinMatches.txt' $peptide_seq_f  $out_dir'proteins2genomesNames_dic.json' $bins2taxonDic_f $out_dir

python3 make_iTOL_node_mostSpecificTaxaAssignment_peptideMapCounts.py  $tree_f $bins2taxonDic_f $allNodes2taxonDic_f $out_dir'peptide2allProteinMatches_allGenomes2Peptides_dic.json' $out_dir'peptide2allProteinMatches_Peptides2allGenomes_dic.json' $out_dir

python3 make_iTOLSymbolFileFromPeptideCountsDic.py $allNodes2taxonDic_f $out_dir'nodeName2NpeptidesMapped_dic.json' $out_dir
