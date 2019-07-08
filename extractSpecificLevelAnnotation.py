# -*- coding: utf-8 -*-
"""
Created on Mon May 27 10:45:17 2019

@author: mstambou

script which takes in the binID 2 the different level of taxa mappings and will 
return a combined dictionary of a specific level binID to taxonomy mappings.
For example return all bin2phylum level mappings.
It will also return another  dictionary mapping that maps between the different phylum 
to different colors, which will later be used in creating a color range file.

this script requires one command line argument which is a dictionary mapping between bin IDs
to different levels of taxonomic assignments, i.e. species, genus, order etc...
second input is the level at which the annotation is to be extractesd i.e. phylum class etc...
third input is the output diretory to store the generated files.
"""

import json
import os
import sys



all_taxa_f = 'bins2taxonomic_assignment_gtdbtk/allBin2taxon_dic_mergedABC_withARCHAEA_dic.json'

if len(sys.argv) != 4:
    print('please enter 3 command line argument1 to run this script, i.e. example to run\n python3 extractSpecifcLevelAnnotation.py bins/2/taxa/map/file annotate_level dir/to/output/directory/')

else:

    all_taxa_f = sys.argv[1]
    #out_dir = hbc_taxa_f.rsplit('/', 1)
    out_dir = sys.argv[3]
    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)
    
    level = sys.argv[2]
    
    with open(all_taxa_f, 'r') as in_f:
        allBin2taxon_dic_mergedABC_withARCHAEA_dic = json.load(in_f)
        
    def getBin2TaxaMap_dic(hbc_taxa_dic, NCBI_taxa_dic, umgs_taxa_dic, level):
        taxa2singleChar_dic = {
                               'kingdom':'d',
                               'phylum':'p',
                               'class':'c',
                               'order':'o',
                               'family':'f',
                               'genus':'g',
                               'species':'s'
                               }
        
        level_char = taxa2singleChar_dic[level]
        
        levelBin2TaxaMap_dic = dict()
        
        for binID in hbc_taxa_dic:
            taxa = hbc_taxa_dic[binID]    
            if level_char in taxa:
                if hbc_taxa_dic[binID][level_char] != '':
                    levelBin2TaxaMap_dic[binID] = hbc_taxa_dic[binID][level_char]
                else:
                    levelBin2TaxaMap_dic[binID] = 'NA'
        
        for binID in NCBI_taxa_dic:
            taxa = NCBI_taxa_dic[binID]
            clean_binID = '_'.join(binID.split('_', 2)[:2])
            if level_char in taxa:
                if NCBI_taxa_dic[binID][level_char] != '':
                    levelBin2TaxaMap_dic[clean_binID ] = NCBI_taxa_dic[binID][level_char]
                else:
                    levelBin2TaxaMap_dic[clean_binID ] = 'NA'
        
        for binID in umgs_taxa_dic:
            taxa = umgs_taxa_dic[binID]
            if level_char in taxa:
                if umgs_taxa_dic[binID][level_char] != '':
                    levelBin2TaxaMap_dic[binID] = umgs_taxa_dic[binID][level_char]
                else:
                    levelBin2TaxaMap_dic[binID] = 'NA'
        
        return levelBin2TaxaMap_dic
    
    def getAllBins2TaxaMap_dic(allBin2taxon_dic_mergedABC_withARCHAEA_dic, level):
        taxa2singleChar_dic = {
                               'kingdom':'d',
                               'phylum':'p',
                               'class':'c',
                               'order':'o',
                               'family':'f',
                               'genus':'g',
                               'species':'s'
                               }
        
        level_char = taxa2singleChar_dic[level]
        
        levelBin2TaxaMap_dic = dict()
        
        for binID in allBin2taxon_dic_mergedABC_withARCHAEA_dic:
            taxa = allBin2taxon_dic_mergedABC_withARCHAEA_dic[binID]    
            if level_char in taxa:
                if allBin2taxon_dic_mergedABC_withARCHAEA_dic[binID][level_char] != '':
                    levelBin2TaxaMap_dic[binID] = allBin2taxon_dic_mergedABC_withARCHAEA_dic[binID][level_char]
                else:
                    levelBin2TaxaMap_dic[binID] = 'NA'    
    
        return levelBin2TaxaMap_dic
    
    
    if __name__ == "__main__":
        
        #levelBin2TaxaMap_dic = getBin2TaxaMap_dic(hbc_taxa_dic, NCBI_taxa_dic, umgs_taxa_dic, level)
        levelBin2TaxaMap_dic = getAllBins2TaxaMap_dic(allBin2taxon_dic_mergedABC_withARCHAEA_dic, level)
        
        all_taxa = list(set(list(levelBin2TaxaMap_dic.values())))
        all_taxa.sort()
        if 'NA' in all_taxa:
            all_taxa.remove('NA')
            all_taxa.append('NA')
        
        color_spectra = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#F09E5F', '#E8B269', '#A5BBFF', '#5C1C21', '#1A2000', '#E8E906', '#92E97C', '#2E0854', '#17644A']                     
        
        if len(all_taxa) > len(color_spectra):
            import random
            color_spectra = list()
            for i in range(len(all_taxa)):
                r = lambda: random.randint(0,255)
                hex_color = ('#%02X%02X%02X' % (r(),r(),r()))
                color_spectra.append(hex_color)
            
            
                         
        taxa2colorMap_dic = dict()
        
        for taxa, color in zip(all_taxa, color_spectra):
            taxa2colorMap_dic[taxa] = color
            
        with open(out_dir+level+'_LevelAllBin2TaxaMap_dic.json', 'w') as out_f:
            json.dump(levelBin2TaxaMap_dic, out_f)
            
        with open(out_dir+level+'_allTaxa2colorMap_dic.json', 'w') as out_f:
            json.dump(taxa2colorMap_dic, out_f)
            
    
        