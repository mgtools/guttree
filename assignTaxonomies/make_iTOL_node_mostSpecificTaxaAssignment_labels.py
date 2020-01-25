# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 18:42:39 2019

@author: mstambou

script in which it will take leaf nodes to taxa assignment dictionary and will create an
iTOL importable file where it assigns all nodes (leaves and internal ones) to the most 
sepcific taxonomical assignment possible
"""

from ete3 import Tree
import json
import sys

tree_in_f = 'allBinsWithArchaea_tree_traditional_annotation_GTDBTK/ribosomal_GTP_EFTU_concatinated_fasttree_internalNodesNamed_pruned_rooted.outtree'
out_dir = tree_in_f.split('/', 1)[0]+'/'
allNodes2taxon_dic_f = 'allBinsWithArchaea_tree_traditional_annotation_GTDBTK/allNodes2taxon_dic.json'

if len(sys.argv) != 4:
    print('Please specify 3 command line arguments, example to run script, python3 make_iTOL_node_mostSPecificTaxaASsignment_labels.py dir/to/tree_f all/nodes/2/taxa/dic_f dir/to/output/')
else:
    tree_in_f = sys.argv[1]
    allNodes2taxon_dic_f = sys.argv[2]
    out_dir = sys.argv[3]
    
    tree = Tree(tree_in_f, format = 1)
    
    with open(allNodes2taxon_dic_f, 'r') as in_f:
        allNodes2taxa_dic = json.load(in_f)
    
    leaf_names = tree.get_leaf_names()    
    all_nodes = list()
    
    for node in tree.iter_prepostorder():
        all_nodes.append(node)
        
    def getNodeName2NodeMap_dic(tree):
        """
        simple function that will return a node name to node mapping
        """
        nodeName2Node_dic = dict()
        for node in tree.iter_prepostorder():
            node = node[1]
            nodeName2Node_dic[node.name] = node
        return nodeName2Node_dic
    
    tree.get_tree_root().name = 'OROOT'
    
    nodeName2Node_dic = getNodeName2NodeMap_dic(tree)
        
    def getMostSpecificTaxon(taxon_dic):
        for char in 'sgfocpd':
            if char in taxon_dic:
                if (taxon_dic[char] != '') and ('GCA' not in taxon_dic[char]) and ('GCF' not in taxon_dic[char]):
                    return char, taxon_dic[char]
        else:
            return 'NA', 'NA'
                
    
    #with open(out_dir+'/'+'allNodes_labels.txt', 'w') as out_f:
    with open(out_dir+'/'+'allNodes_labels_noNumbers.txt', 'w') as out_f, open(out_dir+'/'+'allNodes_labels_withNumbers.txt', 'w') as out_2_f:
        out_f.write('LABELS\n')
        out_2_f.write('LABELS\n')
        out_f.write('SEPARATOR TAB\n')
        out_2_f.write('SEPARATOR TAB\n')        
        out_f.write('\n')
        out_2_f.write('\n')
        out_f.write('DATA\n')
        out_2_f.write('DATA\n')
        for genome in allNodes2taxa_dic:                        
            level, mostSpecificTaxon = getMostSpecificTaxon(allNodes2taxa_dic[genome])
            print(level, mostSpecificTaxon)
            node = nodeName2Node_dic[genome]
            n_leaves = len(node.get_leaves())
            out_2_f.write(genome+'\t'+level+': '+mostSpecificTaxon+' ('+str(n_leaves)+')\n')
            out_f.write(genome+'\t'+mostSpecificTaxon + '\n')
    
    
    with open(out_dir+'/'+'allNodes_popupInfo.txt', 'w') as out_f:
        out_f.write('POPUP_INFO\n')
        out_f.write('\n')
        out_f.write('SEPARATOR TAB\n')
        out_f.write('\n')
        out_f.write('DATA\n')
        for node in allNodes2taxa_dic:
            taxa = allNodes2taxa_dic[node]
            #level, mostSpecificTaxa = getMostSpecificTaxon(taxa)
            popup_text = ''
            popup_text += '<h3 style="font-size:100%;">'+'Bin ID: '+node+'</h3>'
            popup_text += '<h1 style="font-size:100%;">'+'d: '+taxa['d']+'</h1>'
            popup_text += '<h1 style="font-size:100%;">'+'p: '+taxa['p']+'</h1>'
            popup_text += '<h1 style="font-size:100%;">'+'c: '+taxa['c']+'</h1>'
            popup_text += '<h1 style="font-size:100%;">'+'o: '+taxa['o']+'</h1>'
            popup_text += '<h1 style="font-size:100%;">'+'f: '+taxa['f']+'</h1>'
            popup_text += '<h1 style="font-size:100%;">'+'g: '+taxa['g']+'</h1>'
            popup_text += '<h1 style="font-size:100%;">'+'s: '+taxa['s']+'</h1>'
            out_f.write(node+'\tspecies information\t'+popup_text+'\n')
