#!/usr/bin/env Rscript
#R function from the ape package to compute the pairwise distances between the pairs of tips from a phylogenetic tree using it's branch lengths
library(ape)

args = commandArgs(trailingOnly=TRUE)

print(args)

if (length(args) != 2){
   stop("please provide 2 command line arguments to run this script, example to run\n Rscript cophynetic.r dir/to/tree/files/ output/dir/", call.=FALSE)
 }else if (length(args)==2){



full_dir <- file.path(paste(getwd(),'/', args[2], sep = ''))

if (!dir.exists(full_dir))
{
  dir.create(full_dir)
}

tree_files = list.files(args[1])

for (item in tree_files)
{
  in_tree = read.tree(paste(args[1], item, sep = ''))
  pfam = strsplit(item, '_hits_msa')[[1]][1]
  print(pfam)
  leaf_pairwise_dist <- cophenetic(in_tree)
  #pairwise_dist = dist.nodes(in_tree)
  write.table(leaf_pairwise_dist, file =  paste(args[2], pfam, '_pairwiseDistMat.txt', sep = ''), quote = FALSE)
}

}