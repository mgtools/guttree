# -*- coding: utf-8 -*-
"""
Created on Mon May 27 11:58:28 2019

@author: mstambou

script that will get mapping dictionaries between binIDs to taxa, at specified level, mappings
and another dictionary that is taxa level to color mapping, both of which produced by 
'extractSpecificLevelAnnotation.py' script, as a result this script will generate a 
template file that could be dragged and dropped to the iTOL tree to annotate it.

the two command line arguments are the two dictionaries produced by the extractSpecificLevelANnotation.py
third input is the taxonomic level that I want, i.e. phylum, class etc...
fourth input is the output diretory to store the generated files.
"""
import sys
import json

level = 'phylum'
import os

Bin2TaxaMap_f = 'allBinsWithArchaea_tree_averaging_annotation_GTDBTK/'+level+'LevelAllBin2TaxaMap_dic.json'
taxa2colorMap_f = 'allBinsWithArchaea_tree_averaging_annotation_GTDBTK/'+level+'allTaxa2colorMap_dic.json'

if len(sys.argv) != 5:
    print('please enter 4 command line arguments to run this script, i.e. example to run\n python3 make_iTOLcolorStylesFile.py bins/2/taxa/map/file taxa/2/colors/mapfile taxa_level dir/to/output/directory/')

else:
    Bin2_levelTaxaMap_f = sys.argv[1]
    taxa2colorMap_f = sys.argv[2]
    level = sys.argv[3]
    
    #out_dir = Bin2TaxaMap_f.rsplit('/', 1)[0]+ '/'    
    out_dir = sys.argv[4]
    
    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)
    with open(Bin2_levelTaxaMap_f, 'r') as in_f:
        Bin2TaxaMap_dic = json.load(in_f)
        
    with open(taxa2colorMap_f, 'r') as in_f:
        taxa2colorMap_dic = json.load(in_f)
        
        
    ####creating the template file
        
    with open(out_dir+level+'LevelTreeColorRanges.txt', 'w') as out_f:
        out_f.write('TREE_COLORS\n')
        out_f.write('SEPARATOR TAB\n\n')    
        out_f.write('DATA\n')
        for binID in Bin2TaxaMap_dic:
            out_f.write(binID+'\t'+'range'+'\t'+str(taxa2colorMap_dic[Bin2TaxaMap_dic[binID]]+'\n'))
        