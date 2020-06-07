"""
script that will get the fragment 2 genomes table output and phylogenetic tree as input 
and will return for each fragment it's least common ancestory assignment for this read
taking care of unambigious read mapping for the multimapped reads accross the genomes.
the script will also output the number of uniquely mapped reads at each level, i.e. at
each taxonomical level (order class etc...)
"""

from ete3 import Tree
import json
import sys
import pandas as pd
import os
import multiprocessing as mp


sam_f = sys.argv[1]
Bin2taxon_clean_dic_f = sys.argv[2]
allNodes2taxon_dic_f = sys.argv[3]
tree_in_f = sys.argv[4]
n_threads = int(sys.argv[5])
out_dir = sys.argv[6]
CHUNKSIZE = 100000

reads2genomes_f = sam_f.rsplit('.', 1)[0]+'_fragment_out.tsv'

tree = Tree(tree_in_f, format = 1)
tree.get_tree_root().name = 'OROOT'

genome_lst = tree.get_leaf_names()

def filter_frame(df, genome_lst = genome_lst):
    print('processing...')
    for i, row in df.iterrows():
        genomes = [item.replace('.fasta', '') for item in row['genomes'].split('|')]
        genomes = '|'.join(list(set(genomes).intersection(genome_lst)))
        df.loc[i] = [row[0], genomes]
    print('df shape from inside the function', df.shape)
    return df


def getMostSpecificTaxon(taxon_dic):
    taxa_str = list()
    for char in 'dpcofgs':
        if char in taxon_dic:
            if (taxon_dic[char] != '') and ('GCA' not in taxon_dic[char]) and ('GCF' not in taxon_dic[char]):
                taxa_str.append(char+'__'+ taxon_dic[char])
    if len(taxa_str) > 0:
        return '|'.join(taxa_str)
    else:
        return 'unclassified'
def getNodeName2NodeMap_dic(tree):
    """                                                                                                           
    simple function that will return a node name to node mapping                                                  
    """
    nodeName2Node_dic = dict()
    for node in tree.iter_prepostorder():
        node = node[1]
        nodeName2Node_dic[node.name] = node
    return nodeName2Node_dic

nodeName2Node_dic = getNodeName2NodeMap_dic(tree)

#tmp_out = 'tmp_out/'
#if not os.path.isdir(tmp_out):
#    os.mkdir(tmp_out)
def LCA_frame(df, idx, tree = tree, reads2genomes_f = reads2genomes_f):
    print('processing frame', idx)
    LCA_taxa2fragments_dic = dict()
    LCA_nodes2fragments_dic = dict()
    #lca_out_f = open(tmp_out + reads2genomes_f.rsplit('.', 1)[0]+'_LCA_'+str(idx)+'.tsv', 'w')
    #lca_out_f.write('fragment_id\tcommon_ancestor\tLCA_taxa\n')
    for i, row in df.iterrows():
        fragment = row['fragment_id']
        genomes = row['genomes']
        all_ancestors = list()
        if type(genomes) == float:
            continue;
        genomes = [item.replace('.fasta', '') for item in row['genomes'].split('|')]
        if len(genomes) == 1:
            common_ancestor = genomes[0]
            all_ancestors.append(common_ancestor)
            if common_ancestor != 'OROOT':
                node = nodeName2Node_dic[common_ancestor]            
                all_ancestors.extend([item.name for item in  node.get_ancestors()])
        else:
            common_ancestor = tree.get_common_ancestor(genomes).name
            all_ancestors.append(common_ancestor)
            if common_ancestor != 'OROOT':
                node = nodeName2Node_dic[common_ancestor]
                all_ancestors.extend([item.name for item in  node.get_ancestors()])

        if common_ancestor == 'OROOT':
            LCA_taxa = 'd__Bacteria'
        else:
            LCA_taxa = getMostSpecificTaxon(allNodes2taxa_dic[common_ancestor])
        #lca_out_f.write('\t'.join([fragment, common_ancestor, LCA_taxa])+'\n')

        for ancestor in all_ancestors:        
            if ancestor in LCA_nodes2fragments_dic:
                LCA_nodes2fragments_dic[ancestor] += 1
            else:
                LCA_nodes2fragments_dic[ancestor] = 1
    
        for k in range(LCA_taxa.count('|')+1):
            taxa = LCA_taxa.rsplit('|', k)[0]
            if taxa in LCA_taxa2fragments_dic:
                LCA_taxa2fragments_dic[taxa] += 1
            else:
                LCA_taxa2fragments_dic[taxa] = 1
    #lca_out_f.close()
    return LCA_nodes2fragments_dic, LCA_taxa2fragments_dic



if __name__ == '__main__':

    in_df = pd.read_csv(reads2genomes_f, sep = '\t', chunksize = CHUNKSIZE)
    pool = mp.Pool(n_threads)
    print(type(in_df))

    funclist = []

    for df in in_df:
        f = pool.apply_async(filter_frame, [df])
        funclist.append(f)

    out_df = pd.DataFrame()
    all_dfs = list()
    for f in funclist:
        all_dfs.append(f.get())

    out_df = pd.concat(all_dfs)
    out_f = reads2genomes_f.rsplit('.', 1)[0]+'_filtered.tsv'
    print(out_df.shape)
    out_df.dropna(axis = 0, how = 'any', inplace = True)
    print(out_df.shape)
    out_df.to_csv(out_f, sep = '\t', index = False)
    
    reads2genomes_df = pd.read_csv(out_f, sep = '\t', chunksize = CHUNKSIZE)
    
    with open(Bin2taxon_clean_dic_f) as in_f:
        Bin2taxon_clean_dic = json.load(in_f)
    with open(allNodes2taxon_dic_f, 'r') as in_f:
        allNodes2taxa_dic = json.load(in_f)
    
    pool = mp.Pool(n_threads)
    funclist = []
    for i, df in enumerate(reads2genomes_df):
        f = pool.apply_async(LCA_frame, [df, i])
        funclist.append(f)
        
    nodes2fragments_dic_lst, taxa2fragments_dic_lst = list(), list()
    for f in funclist:
        nodes2fragments, taxa2fragments = f.get()
        nodes2fragments_dic_lst.append(nodes2fragments)
        taxa2fragments_dic_lst.append(taxa2fragments)

    
    combined_LCA_nodes2fragments_dic = dict()
    for dic in nodes2fragments_dic_lst:
        for k in dic:
            if k in combined_LCA_nodes2fragments_dic:
                combined_LCA_nodes2fragments_dic[k] += dic[k]
            else:
                combined_LCA_nodes2fragments_dic[k] = dic[k]

    with open(out_dir + sam_f.rsplit('.', 1)[0]+'_LCA_nodes2fragments_dic.json', 'w') as nodes2fragments_dic_out_f:
        json.dump(combined_LCA_nodes2fragments_dic, nodes2fragments_dic_out_f)
    
    combined_LCA_taxa2fragments_dic = dict()
    for dic in taxa2fragments_dic_lst:
        for k in dic:
            if k in combined_LCA_taxa2fragments_dic:
                combined_LCA_taxa2fragments_dic[k] += dic[k]
            else:
                combined_LCA_taxa2fragments_dic[k] = dic[k]
                
    df = pd.DataFrame(columns = ['taxonomy', 'n_fragment'])

    for i, k_v in enumerate(combined_LCA_taxa2fragments_dic.items()):
        df.loc[i] = [k_v[0], k_v[1]]

    df.sort_values(by='taxonomy', inplace =True)

    df.to_csv(out_dir + out_f.rsplit('.', 1)[0]+'_LCA_taxa2fragments.tsv', sep = '\t', index = False, header = None)
