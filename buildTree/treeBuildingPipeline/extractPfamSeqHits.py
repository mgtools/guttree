# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 12:20:41 2019

@author: mstambou

script that goes through the text output resulting from an HMMSCAN and I will define different functions
to extract different information from these output files, such as the start and end of the query to 
target mappings, sequence alignment protions between the target and the query. I can add more functionalties
later in time whenever I need it.

to run this script you require to pass three command line arguments, 
the first one is the directory for the hmmscan out search
the second one is the number of threads to be used by the scripts
the third one is the directory for the output folder where the script will dump it's outputs
"""

from Bio import SearchIO
import os
import sys
import pandas as pd
import shutil
import time
import multiprocessing as mp
import sys

#n_threads = 80
#out_dir = ''
#hmmscan_out_dir = '../genomes.FGS_hmmscan_out/'

if len(sys.argv) != 4:
    print('please specify three command line arguments to run script, example to run the script i.e.\n python3 extractPfamSeqHits.py dir/to/hmmscan_output/ n_threads /dir/to/output/')

else:
    hmmscan_out_dir = sys.argv[1]
    n_threads = int(sys.argv[2])
    out_dir = sys.argv[3]

    repr_genomes_list = [item for item in os.listdir(hmmscan_out_dir) if item.endswith('.txt')]

    print(len(set(repr_genomes_list)))
    print('number of threads used: ' + str(n_threads))

    gene2pfam_hits_dir = 'gene2pfam_hits/'
    gene2bestpfam_hits_dir = 'gene2bestpfam_hits/'

    if not os.path.isdir(out_dir + gene2pfam_hits_dir):    
        os.mkdir(out_dir + gene2pfam_hits_dir) 

    if not os.path.isdir(out_dir + gene2bestpfam_hits_dir):
        os.mkdir(out_dir + gene2bestpfam_hits_dir) 

    def get_binPfamHits(genome_bin):
        """
        basically defnining a function here so that I can parrellize this call using the multiprocessing module in python
        """
        time_start = time.time()
        
        hmm_text_in_f  = hmmscan_out_dir+genome_bin
        
        print('processing genome: '+str(genome_bin))

        with open(hmm_text_in_f, 'r') as in_f:
            hmmscan_hit_df = pd.DataFrame(columns = ['gene_id', 'target_id', 'e_value', 'bitscore', 'gene_hit_seq'])
            for res in SearchIO.parse(in_f, 'hmmer3-text'):
                res_df = pd.DataFrame(columns = ['gene_id', 'target_id', 'e_value', 'bitscore', 'gene_hit_seq'])
                for  hsp in res.hsps:#one gene could have multiple domains, also same profile 
                    hsp_bitscore = hsp.bitscore#could be mapped more than once at different locations to the same gene
                    hsp_eval = hsp.evalue
                    hsp_qseq = str(hsp.query.seq)
                    hsp_qid = hsp.query.id
                    hsp_tid = hsp.hit_id #id of the pfam hit to this particualar gene
                    hsp_qseq = hsp_qseq.replace('-', '').upper()
                    #when reporting the query sequence portion aligned with the pfam profile
                    #remove the gaps, the multiple sequence alignment later will encorporate this.
                    row = [hsp_qid, hsp_tid, hsp_eval, hsp_bitscore, hsp_qseq]
                    hmmscan_hit_df.loc[len(hmmscan_hit_df)] = row
                #res_df = res_df.sort_values(by = ['bitscore'], ascending = False)
                #res_df = res_df.drop_duplicates(subset = ['gene_id', 'target_id'], keep = 'first')
                #frames = [hmmscan_hit_df, res_df]
                #hmmscan_hit_df = pd.concat(frames)
        hmmscan_hit_df = hmmscan_hit_df.sort_values(by = ['bitscore'], ascending = False)
        hmmscan_hit_df = hmmscan_hit_df.drop_duplicates(subset = ['gene_id', 'target_id'], keep = 'first')
        hmmscan_hit_df.to_csv(out_dir + gene2pfam_hits_dir+genome_bin+'##gene2pfam_hits.txt', sep = '\t', index = False)
        best_hmmscan_hit_df = hmmscan_hit_df.sort_values(by = ['bitscore'], ascending = False)
        best_hmmscan_hit_df = best_hmmscan_hit_df.drop_duplicates(subset = ['target_id'], keep = 'first')
        best_hmmscan_hit_df.to_csv(out_dir + gene2bestpfam_hits_dir+genome_bin+'##gene2bestpfam_hits.txt', sep = '\t', index = False)

        print('procesed in '+ str(time.time() - time_start))

    pool = mp.Pool(n_threads)

    zip([*pool.map(get_binPfamHits, repr_genomes_list)]) #I had to put the whole thing ino a list inside the zip since zip expect iterable
    
