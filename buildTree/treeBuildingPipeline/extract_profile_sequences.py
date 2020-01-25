# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 18:38:15 2019

@author: mstambou

script that will take a dictionary that was created from the extract_common_profiles.py
script and will create a multifasta file one for each profile where it contains all the 
sequences comming from the different genes from the different bins for each profile.

takes three command line arguments
second one is the bin to pfam sequences maping dictionary created by the script binID2BestPfamSeqs.py
the third one is the output directory specifying where to store the generated multi fasta files
"""

import json
import os
import shutil
import pandas as pd
import sys

if len(sys.argv) != 3:
    print('please enter 2 command line arguments, example to run i.e. python3 extract_profile_sequences.py /bin2pfam/dic_f dir/to/output/' )

else:
    multi_fasta_dir = sys.argv[2]    

    with open(sys.argv[1], 'r') as in_f:
        pfam2binSeqs_dic = json.load(in_f)
    
    if not os.path.isdir(multi_fasta_dir):
        os.mkdir(multi_fasta_dir)

    for pfam in pfam2binSeqs_dic:
        print(pfam)
        with open(multi_fasta_dir+pfam+'_hits.faa', 'w') as out_f:
            [out_f.write('>'+k+'\n'+v+'\n') for k,v in pfam2binSeqs_dic[pfam].items()]
