# -*- coding: utf-8 -*-    
"""
script that will convert multiple sequence alignments from one format to another
"""

from Bio import AlignIO
import os
import sys

if len(sys.argv) != 3:
    print('please enter 2 command line paramters to run script')

else:
    in_dir = sys.argv[1]
    out_dir = sys.argv[2]
    
    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)
    
    msa_files = os.listdir(in_dir)
    for msa_file in msa_files:
        with open(out_dir + msa_file.rsplit('.', 1)[0]+'.fasta', 'w') as out_f:
            alignment = AlignIO.read(open(in_dir + msa_file), 'stockholm')
            AlignIO.write(alignment, out_f, 'fasta')
