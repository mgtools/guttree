# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 09:35:25 2019

@author: mstambou

simple tree that will take in an unrooted tree and will root it using an outgroup 
"""

from ete3 import Tree

import sys

if len(sys.argv) != 3:
    print('please enter 2 command line arguments, example to run script i.e. python3 rootTree_usingOutgroup.py input/tree/file output/directory/')
else:
    tree_in_f = sys.argv[1]
    out_dir = sys.argv[2]
    
    #archea_genomes = ['GCF_000529525.1', 'GCF_000513315.1', 'GCF_000404225.1', 'GCF_000308215.1', 'GCF_000300255.2', 'GCF_000016525.1', 'GCF_000012545.1']
    #Patescibacteria_genomes = ['UMGS1878', 'UMGS1907',  'UMGS1877', 'UMGS1805', 'UMGS1831', 'UMGS2061', 'UMGS1986', 'UMGS1848', 'UMGS2007']
    outgroup_genomes = list()
    
    genome = input('please enter outgroup genome(bin) IDs without extensions:')
    outgroup_genomes.append(genome)
    while(genome != 'done'):
        outgroup_genomes.append(genome)
        genome = input('please enter outgroup genome(bin) IDs without extensions, type "done" to stop:')
        
        
    def getNodeName2NodeMap_dic(tree):
        """
        simple function that will return a node name to node mapping
        """
        nodeName2Node_dic = dict()
        for node in tree.iter_prepostorder():
            node = node[1]
            nodeName2Node_dic[node.name] = node
        return nodeName2Node_dic         
    
    tree = Tree(tree_in_f, format = 1)   
    
    outgroup_parent_node = tree.get_common_ancestor(outgroup_genomes)
    
    nodeName2Node_dic = getNodeName2NodeMap_dic(tree)
    
    
    tree.set_outgroup(outgroup_parent_node)
    
    root = tree.get_tree_root()
    root.name = 'OROOT'
    
    parent_count = 0
    for node in tree.traverse("postorder"):
        if node.name.startswith('P'):
            count = int(node.name.replace('P', ''))
            if count > parent_count:
                parent_count = count
                
    for node in tree.traverse("postorder"):    
        if node.name == '':
            node.name = 'P'+str(parent_count)
            parent_count += 1
    
    tree.write(outfile = out_dir + tree_in_f.rsplit('/', 1)[-1].rsplit('.', 1)[0]+'_rooted.outtree', format = 1)
