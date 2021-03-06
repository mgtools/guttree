# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 16:54:01 2019

@author: mstambou

script in which it will take leaf nodes to taxa assignment dictionary and will create an
iTOL importable file where it assigns all nodes (leaves and internal ones) to the most 
sepcific taxonomical assignment possible also it will take a dictionary of node name to number 
of peptides mapped to it and will display that in the label as well.

to run this script you need to specify 6 command line arguments:
the first argument is the tree input file in newick format
second argument is the Bin2taxon dictionary mapping file which is produced by 'bins2taxonomic_assignment_GTDBTK.py' script
the third argument is all nodes 2 taxonomic dictionary mapping between all the nodes to taxonomies produced by 'nodes2LCA_maps.py'
the fourth argument is a dictionary mapping between all the genomes 2 unipept peptides dictionary produced by 'peptide2genomeMapping2LCA_stringMatches.py'
the fifth argument is a dictionary mapping between all the peptides to all genomes dictionary produced by 'peptide2genomeMapping2LCA_stringMatches.py'
the sixth argument is the output directory for specifying where to store the files generated by this script

"""

from ete3 import Tree
import json
import sys


#tree_in_f = 'allBins_Archaea_61RibosomalGTP_EFTU_internalNodesNamed_pruned_rooted.outtree'
#Bin2taxon_clean_dic_f = '../bins2taxonomic_assignment_gtdbtk/allBin2taxon_dic_mergedABC_withARCHAEA_dic.json'

#allNodes2taxon_dic_f = '../allBinsWithArchaea_tree_averaging_annotation_GTDBTK/allNodes2taxon_dic.json'

#allGenomes2unipeptS7Peptides_dic_f = 'allGenomes2unipeptS7Peptides_dic.json'
#unipeptS7Peptides2allGenomes_dic_f = 'unipeptS7Peptides2allGenomes_dic.json'

#out_dir = tree_in_f.split('/', 1)[0]+'/'

if len(sys.argv) != 7:
    print('please enter 6 command line arguments to run this script. example to run script i.e. \n python3 make_iTOL_node_mostSpecificTaxaAssignment_peptideMapCounts.py dir/to/tree/file Bin2taxon_clean_dic_f allNodes2taxon_dic_f allGenomes2unipeptS7Peptides_dic_f unipeptS7Peptides2allGenomes_dic_f out_dir')

else:
    tree_in_f = sys.argv[1]
    Bin2taxon_clean_dic_f = sys.argv[2]
    allNodes2taxon_dic_f = sys.argv[3]
    allGenomes2unipeptS7Peptides_dic_f = sys.argv[4]
    unipeptS7Peptides2allGenomes_dic_f = sys.argv[5]
    out_dir = sys.argv[6]
    
    sample_name = 'allBins'
    
    tree = Tree(tree_in_f, format = 1)
    
    with open(Bin2taxon_clean_dic_f) as in_f:
        Bin2taxon_clean_dic = json.load(in_f)
        
    with open(allNodes2taxon_dic_f, 'r') as in_f:
        allNodes2taxa_dic = json.load(in_f)
        
    with open(allGenomes2unipeptS7Peptides_dic_f, 'r') as in_f:
        allGenomes2unipeptS7Peptides_dic = json.load(in_f)
        
    with open(unipeptS7Peptides2allGenomes_dic_f, 'r') as in_f:
        unipeptS7Peptides2allGenomes_dic = json.load(in_f)
        
    total_n_peptides = len(unipeptS7Peptides2allGenomes_dic)
        
    allGenomes2unipeptS7PeptideCounts_dic = {k:len(set(v)) for k,v in allGenomes2unipeptS7Peptides_dic.items()}
    
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
    
    nodeName2NpeptidesMapped_dic = dict()
    leaf_names = tree.get_leaf_names()
    not_found = list()
    
    for nodeName in nodeName2Node_dic:
        if nodeName in leaf_names:
            if nodeName not in allGenomes2unipeptS7Peptides_dic:
                not_found.append(nodeName)
            else:            
                nodeName2NpeptidesMapped_dic[nodeName] = len(set(allGenomes2unipeptS7Peptides_dic[nodeName]))
        else:
            children = nodeName2Node_dic[nodeName].get_leaf_names()
            peptides = list()
            for child in children:
                if child not in allGenomes2unipeptS7Peptides_dic:
                    not_found.append(child)
                else:
                    peptides.extend(allGenomes2unipeptS7Peptides_dic[child])
            peptides = list(set(peptides))
            nodeName2NpeptidesMapped_dic[nodeName] = len(peptides)
        
    def getMostSpecificTaxon(taxon_dic):
        for char in 'sgfocpd':
            if char in taxon_dic:
                if (taxon_dic[char] != '') and ('GCA' not in taxon_dic[char]) and ('GCF' not in taxon_dic[char]):
                    return char, taxon_dic[char]
        else:
            return 'NA', 'NA'
                
    
    with open(out_dir + 'nodeName2NpeptidesMapped_dic.json', 'w') as out_f:
        json.dump(nodeName2NpeptidesMapped_dic, out_f)
    
    #with open(out_dir + sample_name+'_labelsWithoutCounts.txt', 'w') as out_f:
    with open(out_dir + sample_name+'_labels.txt', 'w') as out_f:
        out_f.write('LABELS\n')
        out_f.write('SEPARATOR TAB\n')
        out_f.write('\n')
        out_f.write('DATA\n')
        for i, genome in enumerate(allNodes2taxa_dic):                        
            level, mostSpecificTaxon = getMostSpecificTaxon(allNodes2taxa_dic[genome])
            print(level, mostSpecificTaxon)
            #node = nodeName2Node_dic[genome]
            #n_leaves = len(node.get_leaves())
            n_peptides = 0
            if genome in nodeName2NpeptidesMapped_dic:
                n_peptides = nodeName2NpeptidesMapped_dic[genome]
            out_f.write(genome+'\t'+mostSpecificTaxon+' ('+str(n_peptides)+'/'+str(total_n_peptides)+')\n')
            #out_f.write(genome+'\t'+mostSpecificTaxon+'\n')
    

