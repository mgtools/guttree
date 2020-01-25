# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 17:38:17 2019

@author: mstambou

script that will fill in gaps for the MSAs that miss a particular specie.

"""

from Bio import SeqIO
import os
import sys
import copy

if len(sys.argv) != 4:
    print('please specify 3 command line arguments to run this script, example to run script\npython3 fillMSAgaps.py dir/to/genomes/ MSA/in/dir/ MSA/out/dir/')
else:


    genomes_dir = sys.argv[1]
    MSA_in_dir = sys.argv[2]
    MSA_out_dir = sys.argv[3]

    if os.path.exists(MSA_out_dir) == False:
        os.mkdir(MSA_out_dir)

    all_genome_ids = [item.replace('.fasta', '') for item in os.listdir(genomes_dir)]

    for MSA in os.listdir(MSA_in_dir):
        msa_seqs = list(SeqIO.parse(MSA_in_dir + MSA, 'fasta'))
        msa_len = len(str(msa_seqs[0].seq))
        gaps = '-'*msa_len
        missing_genomes = copy.deepcopy(all_genome_ids)
        with open(MSA_out_dir+MSA.rsplit('.', 1)[0]+'_withGaps.fasta', 'w') as out_f:
            for seq in msa_seqs:
                missing_genomes.remove(seq.id)
                out_f.write('>'+seq.id+'\n')
                out_f.write(str(seq.seq)+'\n')
            for genome in missing_genomes:
                out_f.write('>'+genome+'\n')
                out_f.write(gaps+'\n')
            
    
