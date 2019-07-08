# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 18:35:39 2019

@author: mstambou

script that will create an iTOL importable file that will display symbols on the nodes
whose sizes represent proportions of peptides being mapped at that node level. 

to run this script you require three command line arguments.

the first command line argument is allNodes2taxon_dic dictionary mapping that maps nodes to their
taxonomies that (LCA) which is created by 'nodes2LCA_maps.py'
second command line argument is the node name to number of peptides dictionary mapping 
that is produced by the script 'make_iTOL_node_mostSpecificTaxaAssignment_peptideMapCounts.py'
the third command line input is the output directory for this script to store the generated files.
"""

from ete3 import Tree
import sys
import random
import json
from collections import Counter
import pandas as pd


#out_dir = ''
#sample_name = 'allBins'

if len(sys.argv) != 4:
    print('please specify three command line arguments, example to run this script i.e.\n python3 make_iTOLSymbolFileFromPeptideCountsDic.py dir/2/node2taxon_dic_f dir/2/peptide2genome_dic_f /dir/to/output/file')

else:
    
    allNodes2taxon_dic_f = sys.argv[1]
    nodeName2NpeptidesMapped_dic_f = sys.argv[2]    
    out_dir = sys.argv[3]
    
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
        
    
    with open(nodeName2NpeptidesMapped_dic_f , 'r') as in_f:
        nodeName2NpeptidesMapped_dic = json.load(in_f)
        
    with open(out_dir + sample_name+'_addNodeSymbols.txt', 'w') as out_f:
        out_f.write('DATASET_SYMBOL\n')
        out_f.write('\n')
        out_f.write('SEPARATOR TAB\n ')
        out_f.write('DATASET_LABEL\tpeptide counts labels\n')
        out_f.write('COLOR\t#ffff00\n')
        out_f.write('\n')
        out_f.write('LEGEND_TITLE\tpeptide counts\n')    
        out_f.write('LEGEND_SHAPES\t2\n') #1 for square and 2 for a circle
        out_f.write('LEGEND_COLORS\t#D3D3D3\n')
        out_f.write('LEGEND_LABELS\tmapped peptide\n')        
        out_f.write('\n')
        out_f.write('MAXIMUM_SIZE\t35\n')
        out_f.write('\n')
        out_f.write('DATA\n')
        for genome in nodeName2NpeptidesMapped_dic:
            color = "%06x" % random.randint(0, 0xFFFFFF)
            color = color.upper()
            taxonomicAssignment = allNodes2taxon_dic[genome]
            phylum = getPhylumLevelTaxon(taxonomicAssignment)
            out_f.write(genome+'\t2\t'+str(nodeName2NpeptidesMapped_dic[genome])+'\t#'+color+'\t1\t1\t'+phylum+'\n')
                        
