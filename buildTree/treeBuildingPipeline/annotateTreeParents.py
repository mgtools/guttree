# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 11:36:32 2019

@author: mstambou\

simple script that will go over all the internal parent nodes and will assign names 
to these nodes.

takes two command line arguments, 
first argument is the input tree
second argument is the output directory

"""

from ete3 import Tree
import sys

if len(sys.argv) != 3:
    print('please enter 2 command line arguments, example to run script, i.e.\n python3 annotateTreeParents.py input/tree/file output/directory/')

else:
    tree_in_f = sys.argv[1]
    out_dir = sys.argv[2]
    
    tree = Tree(tree_in_f)
    
    tree.get_leaves()
    
    parent_counter = 1
       
    for node in tree.traverse("postorder"):
        if node.name == '':
            node.name = 'P'+str(parent_counter)
            parent_counter += 1
    
    #tree.name = 'root'                
    tree.write(format = 1, outfile = tree_in_f.split('.')[0]+'_internalNodesNamed.outtree')
    
    open(out_dir+'leaves.txt', 'w').write('\n'.join([item.name for item in tree.get_descendants()]))