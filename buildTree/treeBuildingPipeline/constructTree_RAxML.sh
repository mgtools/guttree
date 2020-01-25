#!/bin/bash
#
#script where I will put things together, all the steps in the pipline will be executed here sequentially by just running this script
#I will have a lot of parameters here however will make them optional for the user if they want to be more specific with them
#
#

START=$(date +%s)
#Default parameter values, changed through command line otherwise
#number of cores used by the programs
n_cores=40
data_dir='../../data/'
tree_dir='../../data/combinedTree/'
distMat_dir='../../data/pfam_FastTree_treeDist/'

while getopts i:t: option
do
case "${option}"
in
    i) genomes_directory=${OPTARG};;
    t) n_cores=${OPTARG};;
esac
done

echo 'number of cores used'
echo $n_cores

echo 'genomes directory'
echo $genomes_directory

genomesFolder=$(echo $genomes_directory | rev | cut -d'/' -f2 | rev)
echo $genomesFolder
echo $data_dir$genomesFolder'.FGS'
FGSFolder=$data_dir$genomesFolder'.FGS/'
HMMSCAN_outFolder=$data_dir$genomesFolder'.FGS_hmmscan_out/'
pfam_db=../treeBuildingData/ribosmal_GTP_EFTU_pfam_db/ribosomal_GTP_EFTU_profiles.hmm


#calling scripts

python3 unifyGenomeExtensions.py $genomes_directory 

sh runFGS_parallel.sh -i $genomes_directory -t $n_cores -o $FGSFolder

sh runHMMSCAN_parallel.sh -i $FGSFolder -t $n_cores -o $HMMSCAN_outFolder -m $pfam_db -e .faa

python3 extractPfamSeqHits.py $HMMSCAN_outFolder $n_cores $data_dir

python3 binID2BestPfamSeqs.py $data_dir'gene2bestpfam_hits/' $data_dir

python3 extract_profile_sequences.py $data_dir'pfam2bins2bestPfamSeqs_dic_of_dics.json' $data_dir'bin2bestPfam_seqs/'

sh get_MSA_parallel.sh -i $data_dir'bin2bestPfam_seqs/' -t $n_cores -o $data_dir'pfam_MSA/'

python3 fillMSAGaps.py $genomes_directory $data_dir'pfam_MSA/' $data_dir'pfam_MSA_withGaps/'

Rscript supermat_phylo.R $data_dir'pfam_MSA_withGaps/' $data_dir'supermat.phylip' $data_dir'gene_partitions.txt'

sed 's/DNA/LG/g' $data_dir'gene_partitions.txt' > $data_dir'protein_partitions.txt'

raxmlHPC -s $data_dir'supermat.phylip' -n 'supermat.tree' -q $data_dir'protein_partitions.txt' -p 12345 -D -T $n_cores -m PROTGAMMAAUTO

mv protein_partitions.txt ../../data/

mv RAxML_* ../../data/

python3 annotateTreeParents.py $data_dir'RAxML_bestTree.supermat.tree' $data_dir
