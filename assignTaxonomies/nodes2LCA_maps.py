# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 09:33:46 2019

@author: mstambou

script that will take all the taxonomic assignments for the bins saved in a dictionary
previously and will propagate that information along the tree and assign least common
ancestory to the internal parent nodes

this script requires two command line arguments to run, the first one is the bins to taxon dictionary 
mapping and the second one being the bacterial tree of life. As as result the script produces a 
dictionary mapping between all the bins (leaves and internal parent nodes) to taxonomical mappings
with the least common ancestors approach.

this script will use two inputs. the first one is a dictionary of bins to taxon mappings produced by the script
'bins2taxonomoc_assignment_GTDBTK.py'
the second input for this script is the phylogenetic gene tree produced in newick format.
third input is the output diretory to store the generated files.
"""

import json
from ete3 import Tree
import copy
import sys

bin2taxon_dic_f = 'bins2taxonomic_assignment_gtdbtk/allBin2taxon_dic_mergedABC_withARCHAEA_dic.json'
tree_in_f = 'allBinsWithArchaea_tree_averaging_annotation_GTDBTK/allBins_Archaea_61RibosomalGTP_EFTU_internalNodesNamed_pruned_rooted.outtree'

if len(sys.argv) != 4:
    print('please enter 3 command line arguments to run this script, i.e. example to run\n python3 nodes2LCA_maps.py bins/2/taxa/map/file dir/to/bacterial/tree/file dir/to/output/directory/')
else:
    bin2taxon_dic_f = sys.argv[1]
    tree_in_f = sys.argv[2]
    out_dir = sys.argv[3]
    
    #out_dir = tree_in_f.rsplit('/')[0]+'/'
    
    with open(bin2taxon_dic_f, 'r') as in_f:
        bin2taxon_dic = json.load(in_f)
        
        
    node2taxon_dic = copy.deepcopy(bin2taxon_dic)
        
    tree = Tree(tree_in_f, format = 1)
    
    leaf_names = tree.get_leaf_names()
    
    node2taxon_dic = {k:v for k,v in node2taxon_dic.items() if k in leaf_names}
    
    visited = tree.get_leaves()
    lookForParents = tree.get_leaves()
    
    node2LCA_dic = dict()
    def checkIfInVisited(child_list, visited):
        """
        simple function that will check if all the children of a certain parent 
        node are allready visited and have a taxonomical assignment
        """
        for child in child_list:
            if child not in visited:
                return False
        return True
    
    def getLCAfromChildren(dic_lst):
        """
        simple function that will take a list of taxonomical assignment for children 
        node and will return the least common ancestor that will explain all the 
        children
        """
        taxa_lst = []
        for taxa in dic_lst:
            taxa_lst.append([k+'__'+v for k,v in taxa.items()])
            LCA_taxa = set.intersection(*map(set, taxa_lst))
        LCA_taxa_dic = {item.split('__')[0]:item.split('__')[-1] for item in LCA_taxa}
        if 'd' not in LCA_taxa_dic:
            LCA_taxa_dic['d'] = ''
        if 'p' not in LCA_taxa_dic:
            LCA_taxa_dic['p'] = ''
        if 'c' not in LCA_taxa_dic:
            LCA_taxa_dic['c'] = ''
        if 'o' not in LCA_taxa_dic:
            LCA_taxa_dic['o'] = ''
        if 'f' not in LCA_taxa_dic:
            LCA_taxa_dic['f'] = ''
        if 'g' not in LCA_taxa_dic:
            LCA_taxa_dic['g'] = ''
        if 's' not in LCA_taxa_dic:
            LCA_taxa_dic['s'] = ''
            
        return LCA_taxa_dic
        
    counter = 0
    while(True):    
        visited_now = []
        print('node2taxon', len(node2taxon_dic))
        print('lookfor parents',len(lookForParents))        
        skipped = []
        for node in lookForParents:
            parent = node.up
            if parent:
                parent_children = parent.children
                if checkIfInVisited(parent_children, visited):
                    if (parent not in visited_now or parent not in visited):
                        parent_name = parent.name
                        children_names = [child.name for child in parent_children]
                        children_taxa = [node2taxon_dic[child] for child in children_names]
                        parent_LCA = getLCAfromChildren(children_taxa)
                        node2taxon_dic[parent_name] = parent_LCA
                        visited_now.append(parent)
                elif node not in skipped:
                    skipped.append(node)
        print('visited now', len(visited_now))
        print('skipped', len(skipped))
        lookForParents = list()
        lookForParents = visited_now
        lookForParents.extend(skipped)
        lookForParents = list(set(lookForParents))
        visited.extend(visited_now)
        if len(lookForParents) == 1 and (tree.get_tree_root() in lookForParents):
            break
        node2taxon_dic['OROOT'] = {'d':'Bacteria', 'p':'', 'c':'', 'o':'', 'f':'', 'g':'', 's':''}
        with open(out_dir + 'allNodes2taxon_dic.json', 'w') as out_f:
            json.dump(node2taxon_dic, out_f)
        