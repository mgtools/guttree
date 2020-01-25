# -*- coding: utf-8 -*-
"""
Created on Tue May  7 16:55:20 2019

@author: mstambou

script that will take as input a matrix dataframe(meaning it only has column names and not rown names)
and will output the same matrix in phylip format, so that I can later use neghborhood joinging to construct 
a gene tree.

this script requires one command line argument to run, which is the combined averaged matrix produced
"""

import pandas as pd
import sys

if len(sys.argv) != 2:
    print('please specify one command line argument to run this script, example to run this script i.e. python3 dfMat2phylip.py dir/to/combined/matrix_f')
else:

    mat_f_in = sys.argv[1]

    mat_df = pd.read_csv(mat_f_in, sep = '\t')
    species = list(mat_df.columns)
    n_species = len(species)
    with open(mat_f_in.rsplit('.', 1)[0]+'.phylip', 'w') as out_f:
        out_f.write(str(n_species)+'\n')
        for i, row in mat_df.iterrows():
            print(species[i])
            str_row = list(map(str, row))
            out_f.write(species[i]+'\t'+'\t'.join(str_row)+'\n')
