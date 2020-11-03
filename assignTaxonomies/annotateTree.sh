#!/bin/bash
#
#script where I will put things together, all the steps needed to annotate a rooted tree starting from a phylogenetic tree in newick format
#and output files of the gtdbtk program
#
#

data_dir='../data/'
out_dir='../assignTaxonomoies/'

taxonomy_level='phylum'


tree_f='../data/combinedTree/allPfamsAveraged_treeDist_clean_internalNodesNamed_rooted.outtree'

while getopts i:o:p:l: option
do
case "${option}"
in
    i) classification_f=${OPTARG};;    
    o) out_dir=${OPTARG};;
    p) tree_f=${OPTARG};;
    l) taxonomy_level=${OPTARG};;
esac
done

data_dir=$out_dir
bins2taxonDic_f=$data_dir'allBin2taxon_dic.json'
allNodes2taxonDic_f=$data_dir'allNodes2taxon_dic.json'

python3 bins2taxonomic_assignment_GTDBTK.py $classification_f $data_dir $tree_f

python3 extractSpecificLevelAnnotation.py $bins2taxonDic_f $taxonomy_level $data_dir

python3 nodes2LCA_maps.py $bins2taxonDic_f $tree_f $data_dir

python3 make_iTOLcolorStylesFile.py $data_dir$taxonomy_level'_LevelAllBin2TaxaMap_dic.json' $data_dir$taxonomy_level'_allTaxa2colorMap_dic.json' $taxonomy_level $out_dir

python3 make_iTOLColorLegendFile.py $data_dir$taxonomy_level'_LevelAllBin2TaxaMap_dic.json' $data_dir$taxonomy_level'_allTaxa2colorMap_dic.json' $taxonomy_level $out_dir

python3 make_iTOL_node_mostSpecificTaxaAssignment_labels.py $tree_f $allNodes2taxonDic_f $out_dir
