# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 17:37:01 2019

@author: mstambou

script to map back the bin names to the newick file created using the kmer distance matrix
since I changed the bin names to numbers with padding of width 10, since otherwise it 
was throwing an error to the phylip package.

this script requires two command line arguments, the first one being the dictionary mapping 
that has mappings from the arbitrary leaf IDs to the actual leaf IDs,
the second one is the newick tree file that has the arbitrary mappings as the leaf IDs.
"""

import json
from Bio import Phylo
import sys

if len(sys.argv) != 3:
    print('please enter two command line arguments to run this script, example to run script i.e.\n python3 mapBackTreeLeafNames.py dir/to/map_dic_f dir/to/newick/tree_f')
else:
    map_dic_f = sys.argv[1]
    tree_f = sys.argv[2]

    with open(map_dic_f, 'r') as in_f:
        map_dic = json.load(in_f)
    
    tree = Phylo.read(tree_f, 'newick')

    for leaf in tree.get_terminals():
        leaf.name = map_dic[leaf.name]
    
    with open(tree_f[:tree_f.rfind('.')]+'_clean.outtree', 'w') as out_f:
        Phylo.write(tree, out_f, 'newick')
