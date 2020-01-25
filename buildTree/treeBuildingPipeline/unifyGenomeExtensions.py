"""
script used to change all the genome extensions to '.fasta'
"""

import os
import sys

if len(sys.argv) != 2:
   print('to run this script please enter 1 command line argument, exampel to run script i.e. python3 unifyGenomeExtensions.py /dir/to/genome/files')

else:
   in_dir = sys.argv[1]
   for file in os.listdir(in_dir):
       out_file = file.rsplit('.', 1)[0]
       os.rename(in_dir+file, in_dir+out_file+'.fasta')

