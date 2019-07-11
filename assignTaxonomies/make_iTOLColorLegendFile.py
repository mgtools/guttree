# -*- coding: utf-8 -*-
"""
Created on Mon May 27 15:10:17 2019

@author: mstambou

script that will also read in the mapping dictionaries produced by the 'extractSpecificLevelAnnotation.py'
script and will create a legend file that could be imported and dragge to the tree

the two command line arguments are the two dictionaries produced by the extractSpecificLevelANnotation.py
third input is the output diretory to store the generated files.
"""

import json
import sys
import os

level = 'phylum'

Bin2TaxaMap_f = 'allBinsWithArchaea_tree_averaging_annotation_GTDBTK/'+level+'LevelAllBin2TaxaMap_dic.json'
taxa2colorMap_f = 'allBinsWithArchaea_tree_averaging_annotation_GTDBTK/'+level+'allTaxa2colorMap_dic.json'


if len(sys.argv) != 5:
    print('please enter 4 command line arguments to run this script, i.e. example to run\n python3 make_iTOLcolorStylesFile.py bins/2/level/taxa/map/file taxa/2/colors/mapfile taxa_level dir/to/output/directory/')

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
        
    taxa2color_pairs = list(taxa2colorMap_dic.items())
    legendShapes = ['1']*len(taxa2color_pairs)    
        
    ####creating the template file
        
    with open(out_dir+level+'LevelTreeLegend.txt', 'w') as out_f:
        out_f.write('DATASET_COLORSTRIP\n')
        out_f.write('SEPARATOR TAB\n\n')   
        out_f.write('DATASET_LABEL\tlabel1\n')
        out_f.write('COLOR\t#FF0000\n\n')
        out_f.write('LEGEND_TITLE\t'+level+' assignments\n')
        out_f.write('LEGEND_SHAPES\t'+'\t'.join(legendShapes)+'\n')
        out_f.write('LEGEND_COLORS\t')
        for item in taxa2color_pairs:
            out_f.write(item[1]+'\t')
        out_f.write('\n')
        out_f.write('LEGEND_LABELS\t')
        for item in taxa2color_pairs:
            out_f.write(item[0]+'\t')
        out_f.write('\n')
        out_f.write('DATA\n')
        for binID in Bin2TaxaMap_dic:
            out_f.write(binID+'\t'+str(taxa2colorMap_dic[Bin2TaxaMap_dic[binID]]+'\t'+Bin2TaxaMap_dic[binID]+'\n'))
