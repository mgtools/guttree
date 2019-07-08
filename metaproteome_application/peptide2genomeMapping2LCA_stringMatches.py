# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 14:20:35 2019

@author: mstambou

script in which I will take exact string matching output and also the peptides fasta file
and also the protein sequence to genome mappings files and will output peptide to genome ID 
dictionary mappings. This dictionary will be later used at another script where it will map 
the peptides to least common ancestor taxonomical assignment.
Also towards the end I will do a least common ancestory mapping of each peptide using the 
allBin2taxon_dic mappings that I calculated over all the genomes using gtdbtk.

the script requires five command line arguments to run
the first one is the output for the peptide matching using string matching generated from the script 'peptide2sequencesStringMatching.py'
the second one is the peptide sequences in fasta format
the third input is a dictionary mapping between protein IDs to their respective genome IDs
the fourth argument is a dictionary mapping of all bins to taxons
the fifth argument is to specify the output directory for thie script to store its outputs
"""

from Bio import SeqIO
import pandas as pd
import json
import sys

#stringMatch_out = 'unipeptS7AllPeptides2myGenomesStringMatch.txt'
#peptide_seqs_f = 'peptideSeqs.fasta'

#proteins2genomes_f = '../test_scripts/data/umgs_hgg_proteins2genomes_dic.json'

#allBin2taxon_f = '../bins2taxonomic_assignment_gtdbtk/allBin2taxon_dic.json'

#out_dir = ''

if len(sys.argv) != 6:
    print('please enter 5 command line arguments, example to run script i.e. \n python3 peptide2genomeMapping2LCA_stringMatches.py dir/2/string/match/file dir/to/peptideseqs_File dir/to/proteins/2/genomes/map_file dir/to/allBin2taxon_dic_file dir/to/output/directory/')

else:
    
    stringMatch_out = sys.argv[1]
    peptide_seqs_f = sys.argv[2]
    proteins2genomes_f = sys.argv[3]
    allBin2taxon_f = sys.argv[4]
    out_dir = sys.argv[5]
    
    stringMatch_out_df = pd.read_csv(stringMatch_out, sep = '\t')
    
    peptide_seqs = SeqIO.to_dict(SeqIO.parse(peptide_seqs_f, 'fasta'))
    
    peptide_seq_lst = [str(peptide_seqs[peptide].seq) for peptide in peptide_seqs]
    
    sample_name = stringMatch_out.split('.')[0]
    
    
    
    with open(proteins2genomes_f, 'r') as in_f:
        proteins2genomesFnames_dic = json.load(in_f)
        
    with open(allBin2taxon_f, 'r') as in_f:
        allBin2taxon_dic = json.load(in_f)
        
    unipeptS7Peptides2allGenomes_dic = dict()
    allGenomes2unipeptS7Peptides_dic = dict()
    
    for i, row in stringMatch_out_df.iterrows():
        peptide_id = row['peptide_ID']
        peptide = str(peptide_seqs[peptide_id].seq)
        #peptide = row['peptide']
        protein_id = row['protein_match_ID']
        genome_id = proteins2genomesFnames_dic[protein_id]
        #if 'GCF_' in genome_id:
        #   genome_id = '_'.join(genome_id.split('_')[:2])
        if peptide in unipeptS7Peptides2allGenomes_dic:
            unipeptS7Peptides2allGenomes_dic[peptide].append(genome_id)
        else:
            unipeptS7Peptides2allGenomes_dic[peptide] = [genome_id]
        if genome_id in allGenomes2unipeptS7Peptides_dic:
            allGenomes2unipeptS7Peptides_dic[genome_id].append(peptide)
        else:
            allGenomes2unipeptS7Peptides_dic[genome_id] = [peptide]
    
    #####map the genomes mapped to peptides to their taxonomic assignments
    unipeptS7Peptides2allGenomes2taxonomies_dic_of_dic = dict()
    unmapped = list()
    
    for peptide in unipeptS7Peptides2allGenomes_dic:
        genomes = unipeptS7Peptides2allGenomes_dic[peptide]
        genomes2taxa_dic = dict()
        for genome in genomes:
            if genome in allBin2taxon_dic:
                taxa = list([str(k)+'__'+str(v) for k,v in allBin2taxon_dic[genome].items()])
                genomes2taxa_dic[genome] = taxa
            else:
                unmapped.append(genome)
        unipeptS7Peptides2allGenomes2taxonomies_dic_of_dic[peptide] = genomes2taxa_dic
    
    unipeptS7Peptides2allGenomes_dic = dict()
    
    for peptide in unipeptS7Peptides2allGenomes2taxonomies_dic_of_dic:
        genome2taxonomies = unipeptS7Peptides2allGenomes2taxonomies_dic_of_dic[peptide]
        genomes = list(genome2taxonomies.keys())
        unipeptS7Peptides2allGenomes_dic[peptide] = list(set(genomes))
        
    with open(out_dir+sample_name+ '_Peptides2allGenomes_dic.json', 'w') as out_f:
        json.dump(unipeptS7Peptides2allGenomes_dic, out_f)
        
    with open(out_dir+sample_name+ '_allGenomes2Peptides_dic.json', 'w') as out_f:
        json.dump(allGenomes2unipeptS7Peptides_dic, out_f)
    
    ###peptide least common ancestory
    unipeptS7Peptides2LCA_dic = dict()
    for peptide in unipeptS7Peptides2allGenomes2taxonomies_dic_of_dic:
        all_taxa = list()
        for genome in unipeptS7Peptides2allGenomes2taxonomies_dic_of_dic[peptide]:
            all_taxa.append(unipeptS7Peptides2allGenomes2taxonomies_dic_of_dic[peptide][genome])
        lca = list(set(all_taxa[0]).intersection(*all_taxa))
        unipeptS7Peptides2LCA_dic[peptide] = lca
        
    
    unipeptS7Peptides2LCA_dic_of_dic = dict()
    
    for peptide in unipeptS7Peptides2LCA_dic:
        taxonomies = unipeptS7Peptides2LCA_dic[peptide]
        taxonomy_dic = {'d':'', 'p':'', 'c':'', 'o':'', 'f':'', 'g':'', 's':''}
        for taxa in taxonomies:
            taxa = taxa.split('__')
            taxonomy_dic[taxa[0]] = taxa[1]
        unipeptS7Peptides2LCA_dic_of_dic[peptide] = taxonomy_dic
        
    with open(out_dir+sample_name+ '_Peptides2LCA_dic_of_dic.json' ,'w') as out_f:
        json.dump(unipeptS7Peptides2LCA_dic_of_dic, out_f)
        
