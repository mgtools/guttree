# -*- coding: utf-8 -*-                                                         
"""                                                                             
Created on Mon Nov 11 17:38:17 2019                                             
                                                                                
@author: mstambou                                                               
                                                                                
script that will convert a multiple sequence alignment file from phylip to fasta
"""     

import os
from Bio import AlignIO
import sys

if len(sys.argv) != 2:
    print('please enter 1 command line arguments tor run this script, example to run script\npython3 phyliptoFasta.py phylip/MSA/file')

else:
    alignment_file = sys.argv[1]
    print(alignment_file)
    with open(alignment_file.rsplit('.', 1)[0]+'.fasta', 'w') as out_f:
        for alignment in AlignIO.read(alignment_file, 'phylip-relaxed'):
            out_f.write('>'+alignment.id+'\n')
            out_f.write(str(alignment.seq)+'\n')
