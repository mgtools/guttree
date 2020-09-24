"""
script that will take a node to LCA nfragments dictionary generated by get_fragment2LCA_node_assignments_report.py
and will create a label file for iTOL annotation
"""

from ete3 import Tree
import json
import sys
import pandas as pd


sam_f = sys.argv[1]
nodes2nfragments_dic_f = sam_f.rsplit('.', 1)[0]+'_LCA_nodes2fragments_dic.json'
Bin2taxon_clean_dic_f = sys.argv[2]
allNodes2taxon_dic_f = sys.argv[3]
out_dir = sys.argv[4]

sample_name = 'allBins'

with open(nodes2nfragments_dic_f) as in_f:
    nodes2nfragments_dic = json.load(in_f)

with open(Bin2taxon_clean_dic_f) as in_f:
    Bin2taxon_clean_dic = json.load(in_f)

with open(allNodes2taxon_dic_f, 'r') as in_f:
    allNodes2taxa_dic = json.load(in_f)

def getMostSpecificTaxon(taxon_dic):
    for char in 'sgfocpd':
        if char in taxon_dic:
            if (taxon_dic[char] != '') and ('GCA' not in taxon_dic[char]) and ('GCF' not in taxon_dic[char]):
                return char, taxon_dic[char]
    else:
        return 'NA', 'NA'



with open(out_dir+sample_name+'_LCA_fragments_labels.txt', 'w') as out_f:
    out_f.write('LABELS\n')
    out_f.write('SEPARATOR TAB\n')
    out_f.write('\n')
    out_f.write('DATA\n')
    for i, node in enumerate(allNodes2taxa_dic):
        level, mostSpecificTaxon = getMostSpecificTaxon(allNodes2taxa_dic[node])
        print(level, mostSpecificTaxon)
        n_fragments = 0
        if node in nodes2nfragments_dic:
            n_fragments = nodes2nfragments_dic[node]
        out_f.write(node +'\t'+mostSpecificTaxon+' ('+str(n_fragments)+')\n')
        