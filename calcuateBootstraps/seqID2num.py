#!/usr/bin/python3

"""
script that will convert the fasta headers into 10 character padded number IDs for phylip interleaved format
"""

from Bio import SeqIO
import sys
import json

def appendZeroes(count):
    s_count = str(count)
    return('0'*(10-len(s_count))+s_count )


if len(sys.argv) != 2:
    print('please enter one command line argument specifying the sequence file in fasta format')
else:
    in_f = sys.argv[1]
    in_seqs = list(SeqIO.parse(in_f, 'fasta'))
    
    out_f = in_f.rsplit('.', 1)[0]+'_paddedIDs.fasta'

    counter2bin_dic = dict()

    with open(out_f, 'w') as output_file:
        for i, seq in enumerate(in_seqs):
            header = seq.id
            sequence = str(seq.seq)
            padded_header = appendZeroes(i)
            counter2bin_dic[padded_header] = header
            output_file.write('>'+padded_header+'\n')
            output_file.write(sequence +'\n')

    with open(in_f.rsplit('.', 1)[0]+'_paddedNum2ID_dic.json', 'w') as out_f:
        json.dump(counter2bin_dic, out_f)
