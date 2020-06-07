"""
script that will output iTOL symbol file
"""

import sys
import random
import json

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

def getPhylumLevelTaxon(taxon_dic):
    if 'p' in taxon_dic:
        if taxon_dic['p'] != '':
            return taxon_dic['p']
    return 'NA'


with open(out_dir + sample_name+'_LCA_fragments_addNodeSymbols.txt', 'w') as out_f:
    out_f.write('DATASET_SYMBOL\n')
    out_f.write('\n')
    out_f.write('SEPARATOR TAB\n ')
    out_f.write('DATASET_LABEL\tLCA fragments labels\n')
    out_f.write('COLOR\t#ffff00\n')
    out_f.write('\n')
    out_f.write('LEGEND_TITLE\tLCA fragments\n')
    out_f.write('LEGEND_SHAPES\t2\n') #1 for square and 2 for a circle                                            
    out_f.write('LEGEND_COLORS\t#D3D3D3\n')
    out_f.write('LEGEND_LABELS\tLCA fragments\n')
    out_f.write('\n')
    out_f.write('MAXIMUM_SIZE\t35\n')
    out_f.write('\n')
    out_f.write('DATA\n')
    for node in allNodes2taxa_dic:
        color = "%06x" % random.randint(0, 0xFFFFFF)
        color = color.upper()
        taxonomicAssignment = allNodes2taxa_dic[node]
        phylum = getPhylumLevelTaxon(taxonomicAssignment)
        n_fragments = 0
        if node in nodes2nfragments_dic:
            n_fragments = nodes2nfragments_dic[node]
        out_f.write(node+'\t2\t'+str(n_fragments)+'\t#'+color+'\t1\t1\t'+phylum+'\n')
