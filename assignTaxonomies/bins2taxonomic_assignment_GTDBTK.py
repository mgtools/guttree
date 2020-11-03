# -*- coding: utf-8 -*-
"""
Created on Sun May 26 14:32:21 2019

@author: mstambou

script that will take the output from gtdbtk classify_wf and will create dictionary 
mappings for all bins taxonomical assignments 
"""

import json
import pandas as pd
import sys
import copy
from ete3 import Tree

if len(sys.argv) != 4:
    print('please specify three command line arguments, example to run script, i.e. python3 bins2taxonomic_assignment_GTDBTK.py bacterial/taxonomic_f archeal/taxonomic_f out/dir')
else:

    taxonomy_f = sys.argv[1]
    out_dir = sys.argv[2]
    tree_f = sys.argv[3]
    
    def get_bin2taxonomyMapping(taxonomy_f, leaves):
        line_count = 0
        bin2taxons_dic = dict()    
            
        taxonomy_df = pd.read_csv(taxonomy_f, sep = '\t')
        for i, row in taxonomy_df.iterrows():
            line_count += 1
            bin_id = row['user_genome']
            if bin_id not in leaves:
                continue
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
    
    tree = Tree(tree_f, format = 1)
    leaves = tree.get_leaf_names()
    allBin2taxon_dic, line_count = get_bin2taxonomyMapping(taxonomy_f, leaves)
    
    kingdom, phylum, _class, order, family, genus, species = get_classificationStats(allBin2taxon_dic)
    
    stats_f = 'allBins_taxaStats.txt'
    dic_f = 'allBin2taxon_dic.json'

    line_count = len(allBin2taxon_dic)
    
    with open(out_dir+stats_f, 'w') as out_f:
        out_f.write('There were '+str(kingdom) + ' kingdom level classification out of '+str(line_count) + ' bins\n')
        out_f.write('There were '+str(phylum) + ' phylum level classification out of '+str(line_count) + ' bins\n')
        out_f.write('There were '+str(_class) + ' class level classification out of '+str(line_count) + ' bins\n')
        out_f.write('There were '+str(order) + ' order level classification out of '+str(line_count) + ' bins\n')
        out_f.write('There were '+str(family) + ' family level classification out of '+str(line_count) + ' bins\n')
        out_f.write('There were '+str(genus) + ' genus level classification out of '+str(line_count) + ' bins\n')
        out_f.write('There were '+str(species) + ' species level classification out of '+str(line_count) + ' bins\n')
    
    for _bin in allBin2taxon_dic:
        taxa = allBin2taxon_dic[_bin]
        if taxa['p'].startswith('Firmicutes'):
            taxa['p'] = 'Firmicutes'
        if taxa['c'].startswith('Bacilli'):
            taxa['c'] = 'Bacilli'
        allBin2taxon_dic[_bin] = taxa
    
    
    with open(out_dir+dic_f, 'w') as out_f:
        json.dump(allBin2taxon_dic,  out_f)
