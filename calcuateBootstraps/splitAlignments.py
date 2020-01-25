#!/usr/bin/python3 

"""
script that will take a phylip formatted multiple sequence alignment file and 
will split it into different files, 1 file per alignment
"""

import sys
from Bio import AlignIO
import os 

if len(sys.argv) != 3:
    print('please enter two command line arguments specifying the phylip formatted alignment file, and the output folder name to create the alignment files')
else:
    in_f = sys.argv[1]
    if not os.path.isdir(sys.argv[2]):
        os.mkdir(sys.argv[2])

    alignments = AlignIO.parse(open(in_f, 'r'), 'phylip')
    
    for i, alignment in enumerate(alignments):
        with open(sys.argv[2]+'alignment_'+str(i+1)+'.phylip', 'w') as out_f:
            AlignIO.write(alignment, out_f, 'phylip')
    
