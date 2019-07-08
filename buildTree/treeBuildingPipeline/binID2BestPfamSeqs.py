# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:07:40 2019

@author: mstambou

script that will output dictionary mappings between each bin to its best pfam sequence hit mappings
that wil later be used in another scritp to extract multi fasta files from them for the MSA.

script takes three command line arguments
the first argument being the list of pfams txt file
second argument is the directory where it contains bin to their best pfam hits
third argugment is the output directory specifying where to dump this mapping
"""

import pandas as pd
import os
import itertools
from operator import itemgetter
import json
import sys

if len(sys.argv) != 3:
    print('please specify 2 comamnd line arguments to run this script, example to run i.e. python3 binID2BestPfamSeqs.py /best/pfam/hits/dir/ out/dir/')

else:          
    best_pfam_hits_dir = sys.argv[1]
    out_dir = sys.argv[2]

    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)

    pfam2bins2bestPfamSeqs_dic_of_dics = dict()
    bin2bestPfam_dic = dict()

    for genome_bin in os.listdir(best_pfam_hits_dir):
        genome_bin_df = pd.read_csv(best_pfam_hits_dir + genome_bin, sep = '\t')
        bin_id = genome_bin.split('.fasta')[0]        
        print(bin_id)
        bin2bestPfam_dic[bin_id] = list(set(genome_bin_df['target_id']))
        for i, row in genome_bin_df.iterrows():
            pfam = row['target_id']
            pfam_seq = row['gene_hit_seq']
            if pfam in pfam2bins2bestPfamSeqs_dic_of_dics:                
                pfam2bins2bestPfamSeqs_dic_of_dics[pfam][bin_id] = pfam_seq
            elif pfam not in pfam2bins2bestPfamSeqs_dic_of_dics:
                pfam2bins2bestPfamSeqs_dic_of_dics[pfam] = dict()
                pfam2bins2bestPfamSeqs_dic_of_dics[pfam][bin_id] = pfam_seq

    with open(out_dir + 'pfam2bins2bestPfamSeqs_dic_of_dics.json', 'w') as out_f:
        json.dump(pfam2bins2bestPfamSeqs_dic_of_dics, out_f)

    bin2Npfam_tup = [(k,len(v) ) for k,v in bin2bestPfam_dic.items()]

    with open(out_dir + 'bin2NPfam.txt', 'w') as out_f:
        out_f.write('bin_id\tNpfams\n')
        for item in bin2Npfam_tup:
            out_f.write(str(item[0])+'\t'+str(item[1])+'\n')
