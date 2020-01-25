#!/usr/bin/env Rscript
#R function from the phylotools package to concatinate the multiple sequence alignments and create a partition file for RAxML in relaxed phylip fomat
library(phylotools)

args = commandArgs(trailingOnly=TRUE)

print(args)

if (length(args) != 3){
   stop("please provide 3 command line arguments to run this script, example to run\n Rscript supermat_phylo.r dir/to/MSAwithGaps/files/ combined/MSA/f_name partition/fname", call.=FALSE)
 }else if (length(args)==3){

full_dir <- file.path(paste(getwd(),'/', args[1], sep = ''))
print(full_dir)
tree_files = list.files(full_dir)
tree_files <- paste(full_dir, tree_files, sep="")
#print(tree_files)
#allowed_tree_files = read.table('../treeBuildingData/allowed_GTDB_120_bac_accessions.txt', sep = '\t', header = FALSE)
#allowed_tree_files <- paste(full_dir, allowed_tree_files$V2, '_hits_msa_withGaps.fasta', sep = "")

supermat(tree_files, outfile = args[2], partition.file = args[3])
#supermat(allowed_tree_files, outfile = args[2], partition.file = args[3])
#supermat(tree_files, outfile = test_supermat.phylip, partition.file = test_partition.txt)

}