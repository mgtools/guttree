# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 18:35:39 2019

@author: mstambou

script that will create iTOL symbol file reflecting the abundance values from metagenomic samples
This script requires two command line arguments, the first one is a dictionary mapping between all 
nodes (both leasves and internal nodes) made from the script 'nodes2LCA_maps.py' and the second 
argument is a dictionary mapping between the name of the nodes (both leaves and internal parent nodes)
to relative abundance as as percentage values dictionary that is produced by the script 
'make_iTOL_node_mostSpecificTaxaAssignment_MGAbundance.py'
"""

import sys
import random
import json

allNodes2taxon_dic_f = '../allBinsWithArchaea_tree_averaging_annotation_GTDBTK/allNodes2taxon_dic.json'
nodeName2MGPercentageAbundance_dic_f = 'nodeName2MGPercentageAbundance_dic.json'

if len(sys.argv) != 4:
    print('please enter 3 command line arguments to run this script, i.e. example to run\n python3 make_iTOLSymbolFileFromMGAbundanceDic.py dir/to/allNodes/toTaxon/file dir/to/nodeName2MGPercentageAbundance_dic_f dir/to/output/directory/')

else:
    out_dir = sys.argv[3]
    allNodes2taxon_dic_f = sys.argv[1]
    nodeName2MGPercentageAbundance_dic_f = sys.argv[2]
    sample_name = 'allBins'
    
    with open(allNodes2taxon_dic_f, 'r') as in_f:        
        allNodes2taxon_dic = json.load(in_f)
        
    def getMostSpecificTaxon(taxon_dic):
        for char in 'sgfocpd':
            if char in taxon_dic:
                if (taxon_dic[char] != '') and ('GCA' not in taxon_dic[char]) and ('GCF' not in taxon_dic[char]):
                    return taxon_dic[char]
            
            
    def getPhylumLevelTaxon(taxon_dic):
        if 'p' in taxon_dic:
            if taxon_dic['p'] != '':
                return taxon_dic['p']
        return 'NA'    
    
    with open(nodeName2MGPercentageAbundance_dic_f, 'r') as in_f:
        nodeName2MGPercentageAbundance_dic = json.load(in_f)
    
    with open(sample_name+'_addNodeSymbols.txt', 'w') as out_f:
        out_f.write('DATASET_SYMBOL\n')
        out_f.write('\n')
        out_f.write('SEPARATOR TAB\n ')
        out_f.write('DATASET_LABEL\tMG abundance labels\n')
        out_f.write('COLOR\t#ffff00\n')
        out_f.write('\n')
        out_f.write('LEGEND_TITLE\tMG abundance\n')    
        out_f.write('LEGEND_SHAPES\t2\n') #1 for square and 2 for a circle
        out_f.write('LEGEND_COLORS\t#D3D3D3\n')
        out_f.write('LEGEND_LABELS\tMG abundance\n')        
        out_f.write('\n')
        out_f.write('MAXIMUM_SIZE\t35\n')
        out_f.write('\n')
        out_f.write('DATA\n')
        for genome in nodeName2MGPercentageAbundance_dic:
            color = "%06x" % random.randint(0, 0xFFFFFF)
            color = color.upper()
            taxonomicAssignment = allNodes2taxon_dic[genome]
            phylum = getPhylumLevelTaxon(taxonomicAssignment)
            out_f.write(genome+'\t2\t'+str(nodeName2MGPercentageAbundance_dic[genome])+'\t#'+color+'\t1\t1\t'+phylum+'\n')
                        
