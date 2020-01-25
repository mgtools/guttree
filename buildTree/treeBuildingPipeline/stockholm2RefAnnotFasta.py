#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 13:00:55 2019

@author: moses

script that will read a stockholm formated multiple sequence alignment resulting 
from hmmalign and will output the MSA in Fasta Format by keeping the parts that 
only align with match states in the HMM model. Guided by the refernce annotation
"""

from Bio import AlignIO
import os
import sys

if len(sys.argv) != 3:
    print('please specify 2 command line arguments to run this script')
else:
    alignment_dir = sys.argv[1]
    out_dir = sys.argv[2]    
    
    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)   
    
    for seq_f in os.listdir(alignment_dir):
        
        alignment = AlignIO.read(alignment_dir + seq_f, 'stockholm')
        
        ref_annot = (alignment.column_annotations)['reference_annotation']
        
        indices = [i for i, char in enumerate(ref_annot) if char == 'x']
        
        with open(out_dir + seq_f.rsplit('.', 1)[0]+'_profileAln.fasta', 'w') as out_f:
            for seq in alignment:
                header = seq.id
                aln = str(seq.seq)
                profile_aln = ''.join([aln[idx] for idx in indices])
                out_f.write('>'+header+'\n')
                out_f.write(profile_aln+'\n')
