#!/usr/bin/python3

from Bio import AlignIO
import sys

if len(sys.argv) != 2:
    print('please enter one comman line argument specifying the fasta formatted sequence file')
else:
    in_f = sys.argv[1]
    
    alignments = AlignIO.read(open(in_f, 'r'), 'fasta')
    AlignIO.write(alignments, open(in_f.rsplit('.', 1)[0]+'.phylip', 'w'), 'phylip')

with open('seqboot_cmds.txt', 'w') as out_f:
    out_f.write(in_f.rsplit('.', 1)[0]+'.phylip\n')
    out_f.write('Y\n')
    out_f.write('173\n')
