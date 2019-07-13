# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:04:33 2019

@author: mstambou

script that will map all protein IDs to genome file names for all the binned and references genomes ~3000 bins
"""

import json
from Bio import SeqIO
import os
import sys

#proteins_dir = 'all_proteins/'
#extension = '.fa.FGS.faa'
#out_dir = 'data/'


if len(sys.argv) != 4:
    print('please enter 3 command line arguments, example to run i.e. python3 proteins2genomes.py dir/2/proteins/ file_extensions dir/2/output/directory/')

else:
    
    proteins_dir = sys.argv[1]
    extension = sys.argv[2]
    out_dir = sys.argv[3]
    proteins2genomes_dic = dict()
    proteins2genomesWithPrefix_dic = dict()
    uncommon = 0
    
    for genome in [file for file in os.listdir(proteins_dir) if file.endswith(extension)]:
        print(genome)
        seqs = list(SeqIO.parse(proteins_dir+genome, 'fasta'))
        genome_name = genome.split(extension)[0]
        for seq in seqs:
            if seq.id in proteins2genomes_dic:
                print(seq.id + ' not a commom sequence ID across genomes')
                uncommon += 1
            proteins2genomes_dic[seq.id] = genome_name
            proteins2genomesWithPrefix_dic[genome_name+':'+seq.id] = genome_name
    
    with open(out_dir + 'proteins2genomesNames_dic.json', 'w') as out_f:
        json.dump(proteins2genomes_dic, out_f)
    
    with open(out_dir + 'proteins2genomesNamesWithPrefix_dic.json', 'w') as out_f:
        json.dump(proteins2genomesWithPrefix_dic, out_f)
