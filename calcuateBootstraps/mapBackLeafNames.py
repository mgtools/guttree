#!/usr/bin/python3

"""
script that will map back the names of the leaves using the dictionary created by the seqID2num.py script.


"""

import sys
from ete3 import Tree
import json

if len(sys.argv) != 3:
    print('please enter 2 command line arguments to run this script, tree file in newick format, dictionary map file')

else:
    tree_f = sys.argv[1]
    map_f = sys.argv[2]

    with open(tree_f, 'r') as in_f1:
        lines = [item.strip() for item in in_f1.readlines()]

    with open(map_f, 'r') as in_f2:
        map_dic = json.load(in_f2)

    tree_out_f = tree_f.replace('paddedIDs', 'originalIDs')

    with open(tree_out_f, 'w') as out_f:
        for line in lines:
            t = Tree(line)
            for node in t.traverse('postorder'):
                if node.is_leaf():
                    node.name = map_dic[node.name]

            out_f.write(t.write()+'\n')
