# -*- coding: utf-8 -*-
"""
Created on Wed May  1 10:18:36 2019

@author: mstambou

script in python that goes over all the individual matrices and then creates one
big matrix that contains all the bins/species from the different individual pfam 
matrices, it does that by averaging out the values for pairs that occur in more
than one pfam matrix.

this script requries one command line argument, the directory for the pfam matrices
"""

import os
from collections import Counter
import pandas as pd
import numpy as np
import sys

if len(sys.argv) != 2:
    print('please specify one command line argument, example to run script python3 combTreeDistMats.py dir/to/pfam/mats/')

else:
    mat_dir = sys.argv[1]

    all_genomes = []

    for mat in os.listdir(mat_dir):
        with open(mat_dir+mat, 'r') as in_f:
            first_line = (in_f.readline()).strip('\n')
            first_line = first_line.split(' ')        
            all_genomes.extend(first_line)
        
    all_genomes_counts = Counter(all_genomes) #this step sees how many matrices each genome is rpesent in

    unique_genomes = list(set(all_genomes)) #list of the set of genomes

    genome2index_dic = {genome:idx for idx, genome in enumerate(unique_genomes)}
    index2genome_dic = {idx:genome for genome,idx in genome2index_dic.items()}

    averagedMSA_pairwiseDistMat_df = pd.DataFrame(columns = unique_genomes, index = unique_genomes)

    pairwiseDist_dic = {genome:dict() for genome in unique_genomes}

    for genome1 in unique_genomes:
        for genome2 in unique_genomes:
            pairwiseDist_dic[genome1][genome2] = []        #creates and fills a dictionary of pairwise distances

    for treeDist_mat_f in os.listdir(mat_dir):
        print('processing: ', treeDist_mat_f, ' ...')
        treeDist_mat_df = pd.read_csv(mat_dir + treeDist_mat_f, sep = ' ')
        treeDist_mat = treeDist_mat_df.as_matrix(columns = list(treeDist_mat_df.columns))
        treeDist_bin2idx_dic = {col:i for i, col in enumerate(treeDist_mat_df.columns)}
        bins = list(treeDist_bin2idx_dic.keys())
        for i, bin1 in enumerate(bins):
            for bin2 in bins[i:]:
                dist_val = treeDist_mat[treeDist_bin2idx_dic[bin1]][treeDist_bin2idx_dic[bin2]]
                pairwiseDist_dic[bin1][bin2].append(dist_val)
                pairwiseDist_dic[bin2][bin1].append(dist_val) #since they are symmetric matrices.
    
    print('calculating final average distance matrix...')
    for i, bin1 in enumerate(unique_genomes):
        print(i, bin1)
        for bin2 in unique_genomes[i:]:
            row, col = genome2index_dic[bin1], genome2index_dic[bin2]
            avg_dist = np.mean(pairwiseDist_dic[bin1][bin2])
            averagedMSA_pairwiseDistMat_df[bin1][bin2] = avg_dist
            averagedMSA_pairwiseDistMat_df[bin2][bin1] = avg_dist
    averagedMSA_pairwiseDistMat_df.to_csv(mat_dir+'allPfamsAveraged_treeDist.txt', sep = '\t', index = False)
