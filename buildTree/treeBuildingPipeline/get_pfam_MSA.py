# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 18:50:21 2019

@author: mstambou

script that given an input for a pfam file will construct a multiple sequence allignment
for that pfams sequences.
"""

from Bio.Align.Applications import MuscleCommandline
import time
import os
import sys

pfam_seq_f = sys.argv[1]
msa_out_dir = sys.argv[2]

#multi_fasta_dir = 
#msa_out_dir = 'allBins_genTree_traditional/MSA_out/'

if os.path.exists(msa_out_dir) == False:
    os.mkdir(msa_out_dir)

def get_MSA(pfam_seq_f):
    msa_out_f = (pfam_seq_f.split('/')[-1]).split('.')[0]+'_msa.fasta'
    #msa_out_f = (pfam_seq_f.split('/')[-1]).split('.')[0]+'_msa.clw'
    
    print('performing multiple sequence alignment')
    start_time = time.time()
    muscle_exe = r"muscle3.8.31_i86linux64"
    muscle_cline = MuscleCommandline(muscle_exe, input=pfam_seq_f, out = msa_out_dir + msa_out_f)
    #muscle_cline = MuscleCommandline(muscle_exe, input=pfam_seq_f, out = msa_out_dir + msa_out_f, clw = True)    
    print(muscle_cline)
    assert os.path.isfile(muscle_exe), 'muscle executable missing'
    stdout, stderr = muscle_cline()
    print('time taken to build MSA '+ str(time.time() - start_time))
    
get_MSA(pfam_seq_f)
