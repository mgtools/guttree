# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 12:13:32 2019

@author: mstambou

script that will take a multiple sequence alignment in some format and will use 
FastTree phylogenetic tree constructing software to build a gene tree.
"""

from Bio.Phylo.Applications import _Fasttree
import time
import os
import sys

print(sys.argv)

msa_seq_f = sys.argv[1]

msa_dir = '/'.join(msa_seq_f.split('/')[:-1])
#tree_out_dir = 'pfamGeneSeqsMSA_Tree/'
tree_out_dir = sys.argv[2]

if os.path.exists(tree_out_dir) == False:
    os.mkdir(tree_out_dir)

def get_geneTree(msa_seq_f):
    
    out_tree = (msa_seq_f.split('/')[-1]).split('.')[0]+'_fasttree.tree'
    print('constructing a phylogenetic tree using FastTree')
    print(msa_seq_f)
    print(tree_out_dir+out_tree)
    start_time = time.time()
    #fasttree_exe = r'FastTree'
    #fasttree_cline = _Fasttree.FastTreeCommandline(fasttree_exe, input=msa_seq_f, out=tree_out_dir+out_tree)
    #print(fasttree_cline)
    #assert os.path.isfile(fasttree_exe), "fasttree executable missing"
    #stdout, stderr = fasttree_cline()
    os.system("FastTree -wag -gamma -pseudo -spr 4 -mlacc 2 -slownni " + msa_seq_f + " > " + tree_out_dir + out_tree)
    print('Tree built in: '+str(time.time() - start_time))

get_geneTree(msa_seq_f)
