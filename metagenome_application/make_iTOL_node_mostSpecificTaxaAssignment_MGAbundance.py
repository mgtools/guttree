# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 16:54:01 2019

@author: mstambou

script in which it will take leaf nodes to taxa assignment dictionary and will create an
iTOL importable file where it assigns all nodes (leaves and internal ones) to the most 
sepcific taxonomical assignment possible also it will take a dictionary of genome to number 
of reads mapped to it and will display that in the label as well.

this script requires five command line arguments.
the first argument is the bacterial tree of life in newick format
the second one is the dictionary mapping between all the bins to their taxonomical assignments
the third argument is a dictionary mapping between all nodes (leaves and internal parent nodes)+
to their least common ancestors, this dictionary is produced by 'nodes2LCA.py' script
and the fourth argument is the abundance table file which is produced by the script called 
'samRead2multiFastaGenomeABundance.py' script
and the fifth argument is the output directory to store the generated files.

This script will generate a dictionary that is a mapping between the node names and the relative 
abundances of each node in the sample as a percentage.
Also it will generate an iTOL importable label file that could be dragged and dropped to the iTOL tree
which will display the percentage abundances of each of the nodes in the sample
"""

from ete3 import Tree
import json
import sys
import pandas as pd

#tree_in_f = 'allBins_Archaea_61RibosomalGTP_EFTU_internalNodesNamed_pruned_rooted.outtree'
#Bin2taxon_clean_dic_f = '../bins2taxonomic_assignment_gtdbtk/allBin2taxon_dic_mergedABC_withARCHAEA_dic.json'

#allNodes2taxon_dic_f = '../allBinsWithArchaea_tree_averaging_annotation_GTDBTK/allNodes2taxon_dic.json'

#abundance_f = 'SRR769523_abundances.txt'

if len(sys.argv) != 6:
    print('please enter 5 command line arguments, to run this script, i.e. example to run\n python3 make_iTOL_node_mostSpecificTaxaASsignment_MGAbundance.py dir/to/tree_f dir/to/bin2taxon/dic_f dir/to/allNodes2taxon/dic_f dir/to/abundanceTable_f dir/2/output/directory/')

else:
    sample_name = 'allBins'
    
    tree_in_f = sys.argv[1]
    Bin2taxon_clean_dic_f = sys.argv[2]
    allNodes2taxon_dic_f = sys.argv[3]
    abundance_f = sys.argv[4]
    out_dir = sys.argv[5]
    
    
    
    abundance_df = pd.read_csv(abundance_f, sep = '\t')
    
    tree = Tree(tree_in_f, format = 1)
    
    with open(Bin2taxon_clean_dic_f) as in_f:
        Bin2taxon_clean_dic = json.load(in_f)
        
    with open(allNodes2taxon_dic_f, 'r') as in_f:
        allNodes2taxa_dic = json.load(in_f)
        
    
    
    leaf_names = tree.get_leaf_names()    
    all_nodes = list()
    
    for node in tree.iter_prepostorder():
        all_nodes.append(node)
        
    def getNodeName2NodeMap_dic(tree):
        """
        simple function that will return a node name to node mapping
        """
        nodeName2Node_dic = dict()
        for node in tree.iter_prepostorder():
            node = node[1]
            nodeName2Node_dic[node.name] = node
        return nodeName2Node_dic
    
    tree.get_tree_root().name = 'OROOT'
    
    nodeName2Node_dic = getNodeName2NodeMap_dic(tree)
     
    def getMostSpecificTaxon(taxon_dic):
        for char in 'sgfocpd':
            if char in taxon_dic:
                if (taxon_dic[char] != '') and ('GCA' not in taxon_dic[char]) and ('GCF' not in taxon_dic[char]):
                    return char, taxon_dic[char]
        else:
            return 'NA', 'NA'
                
    bin2abundance_dic = dict()
    
    for i, row in abundance_df.iterrows():
        genome = row['genome_id']
        total_abundance = row['total_abundance_pmkb']
        genome = genome.rsplit('.', 1)[0]
        bin2abundance_dic[genome] = total_abundance
    
    total_abundance = sum(bin2abundance_dic.values())
    
    bin2percentageAbundance_dic = {k:(v*100)/total_abundance for k,v in bin2abundance_dic.items()}
    
    not_found = list()
    nodeName2MGPercentageAbundance_dic = dict()
    
    for nodeName in nodeName2Node_dic:
        if nodeName in leaf_names:
            if nodeName not in bin2percentageAbundance_dic:
                not_found.append(nodeName)
            else:
                nodeName2MGPercentageAbundance_dic[nodeName] = bin2percentageAbundance_dic[nodeName]
        else:
            children = nodeName2Node_dic[nodeName].get_leaf_names()
            abundances = list()
            for child in children:
                if child not in bin2percentageAbundance_dic:
                    not_found.append(child)
                else:
                    abundances.append(bin2percentageAbundance_dic[child])
            
            nodeName2MGPercentageAbundance_dic[nodeName] = sum(abundances)
    
    
    with open(out_dir+'nodeName2MGPercentageAbundance_dic.json', 'w') as out_f:
        json.dump(nodeName2MGPercentageAbundance_dic, out_f)
    
    #with open(sample_name+'_labelsWithoutCounts.txt', 'w') as out_f:
    with open(out_dir+sample_name+'_labels.txt', 'w') as out_f:
        out_f.write('LABELS\n')
        out_f.write('SEPARATOR TAB\n')
        out_f.write('\n')
        out_f.write('DATA\n')
        for i, genome in enumerate(allNodes2taxa_dic):                        
            level, mostSpecificTaxon = getMostSpecificTaxon(allNodes2taxa_dic[genome])
            print(level, mostSpecificTaxon)
            #node = nodeName2Node_dic[genome]
            #n_leaves = len(node.get_leaves())
            MGPercetageAbundance = 0
            if genome in nodeName2MGPercentageAbundance_dic:
                MGPercetageAbundance = nodeName2MGPercentageAbundance_dic[genome]
            out_f.write(genome+'\t'+mostSpecificTaxon+' ('+str(round(MGPercetageAbundance,3))+'%)\n')
            #out_f.write(genome+'\t'+mostSpecificTaxon+'\n')
    
    
