# -*- coding: utf-8 -*-
"""
Created on Sun May 26 14:32:21 2019

@author: mstambou

script that will take the output from gtdbtk classify_wf and will create dictionary 
mappings for all bins taxonomical assignments 

to run this script you only need to pass the gtdbtk summary.tsv output file.
this will create a dictionary of all the bins to their taxonomy mappings 
which will be used as input to the nodes2LCA_maps.py. Note if you obtain many
GTDBTK files from running multiple instances of GTDBTK then you have to run this 
script on the different GRDBTK output files produced and combine the dictionaries.

also the second command line argument is the directory where you want to store the output
"""

import json
import pandas as pd
import sys

taxonomy_f = 'bins2taxonomic_assignment_gtdbtk/ARCHAEA.ar122.summary.tsv'

if len(sys.argv) != 3:
    print('please enter 2 command line argument1 to run this script, i.e. example to run\n python3 bin2taxonomic_assignment_GTDBTK.py dir/to/GTDBK/out/file.summary.tsv dir/to/output/directory/')

else:
    
    def get_bin2taxonomyMapping(taxonomy_f):
        line_count = 0
        bin2taxons_dic = dict()    
            
        taxonomy_df = pd.read_csv(taxonomy_f, sep = '\t')
        for i, row in taxonomy_df.iterrows():
            line_count += 1
            bin_id = row['user_genome']
            taxonomies = row['classification'].split(';')
            taxonomy_hierarchy = {taxonomy.split('__')[0]:taxonomy.split('__')[-1] for taxonomy in taxonomies}
            bin2taxons_dic[bin_id] = taxonomy_hierarchy
                
        return bin2taxons_dic, line_count
    
    def get_classificationStats(bin2taxons_dic):
        """
        simple function to report how much of the bins were classified to what levels
        """
        kingdom, phylum, _class, order, family, genus, species = 0,0,0,0,0,0,0
        for _bin in bin2taxons_dic:        
            if bin2taxons_dic[_bin]['d'] != '':
                kingdom += 1
            if bin2taxons_dic[_bin]['p'] != '':
                phylum += 1
            if bin2taxons_dic[_bin]['c'] != '':
                _class += 1
            if bin2taxons_dic[_bin]['o'] != '':
                order += 1
            if bin2taxons_dic[_bin]['f'] != '':
                family += 1
            if bin2taxons_dic[_bin]['g'] != '':
                genus += 1
            if bin2taxons_dic[_bin]['s'] != '':
                species += 1
        return kingdom, phylum, _class, order, family, genus, species
    
    
    bin2taxons_dic, line_count = get_bin2taxonomyMapping(taxonomy_f)
    kingdom, phylum, _class, order, family, genus, species = get_classificationStats(bin2taxons_dic)
    stats_f = taxonomy_f.rsplit('/', 1)[-1].split('.')[0]+'_taxaStats.txt'
    dic_f = taxonomy_f.rsplit('/', 1)[-1].split('.')[0]+'Bin2taxon_dic.json'
    #out_dir = taxonomy_f.rsplit('/', 1)[0]+'/'
    out_dir = sys.argv[2]
    with open(out_dir+stats_f, 'w') as out_f:
        out_f.write('There were '+str(kingdom) + ' kingdom level classification out of '+str(line_count) + ' bins\n')
        out_f.write('There were '+str(phylum) + ' phylum level classification out of '+str(line_count) + ' bins\n')
        out_f.write('There were '+str(_class) + ' class level classification out of '+str(line_count) + ' bins\n')
        out_f.write('There were '+str(order) + ' order level classification out of '+str(line_count) + ' bins\n')
        out_f.write('There were '+str(family) + ' family level classification out of '+str(line_count) + ' bins\n')
        out_f.write('There were '+str(genus) + ' genus level classification out of '+str(line_count) + ' bins\n')
        out_f.write('There were '+str(species) + ' species level classification out of '+str(line_count) + ' bins\n')
    
    with open(out_dir+dic_f, 'w') as out_f:
        json.dump(bin2taxons_dic,  out_f)