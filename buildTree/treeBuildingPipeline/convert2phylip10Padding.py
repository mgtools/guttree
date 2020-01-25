# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 14:46:03 2019

@author: mstambou

phylip tree building software for some reason is not able to read the phylip format matrices
that have column IDs that are bigger than 10 for that reason I'm making this script
that will add arbitrary 10 padded number IDs to the species IDs and maintain a dictionary
mapping of this so later it can be reversed in the newick tree file.
"""

import json
import sys

if len(sys.argv) != 3:
    print('please enter two command line argument specifying the distance matrix in phylip format, and the output directory to store the dictionary mapping, example to run the script, python3 convert2phylip10Padding.py /dir/to/distMatf /out/dir/')
else:
    in_mat = sys.argv[1]
    out_mat = (in_mat.rsplit('.', 1)[0] + '_padded.phylip').rsplit('/', 1)[-1]
    out_dir = sys.argv[2]

    def appendZeroes(count):
        s_count = str(count)
        return('0'*(10-len(s_count))+s_count )
    
    counter2bin_dic = dict()

    with open(in_mat, 'r') as in_f, open(out_mat, 'w') as out_f:
        for i, line in enumerate(in_f):
            line = line.split('\t')
            if len(line) == 1:
                new_l = '\t'.join(line)
            else:
                bin_name = line[0]
                bin_number = appendZeroes(i)
                counter2bin_dic[bin_number] = bin_name#.replace(suffix, '')
                new_l = '\t'.join(line[1:])
                new_l = bin_number+'\t'+new_l
            out_f.write(new_l.strip('\n')+'\n')

    with open(out_dir + 'allPfamsAveraged_treeDist_padded_number2bin_dic.json', 'w') as out_f:
        json.dump(counter2bin_dic, out_f)
    
    with open('neighbor_cmds.txt', 'w') as out_f:
        out_f.write('allPfamsAveraged_treeDist_padded.phylip\n')
        out_f.write('Y\n')
