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

sh runFGS.sh -i $genomes_directory -t $n_cores -o $FGSFolder

sh runHMMSCAN_parallel.sh -i $FGSFolder -t $n_cores -o $HMMSCAN_outFolder -m $pfam_db -e .faa

python3 extractPfamSeqHits.py $HMMSCAN_outFolder $n_cores $data_dir

python3 binID2BestPfamSeqs.py $data_dir'gene2bestpfam_hits/' $data_dir

python3 extract_profile_sequences.py $data_dir'pfam2bins2bestPfamSeqs_dic_of_dics.json' $data_dir'bin2bestPfam_seqs/'

sh get_MSA_parallel.sh -i $data_dir'bin2bestPfam_seqs/' -t $n_cores -o $data_dir'pfam_MSA/'

sh get_FastTree_parallel.sh -i $data_dir'pfam_MSA/' -t $n_cores -o $data_dir'pfam_FastTree/'

Rscript cophenetic_phylo.R $data_dir'pfam_FastTree/' $data_dir'pfam_FastTree_treeDist/'

python3 combTreeDistMats.py $data_dir'pfam_FastTree_treeDist/'

python3 dfMat2phylip.py $data_dir'pfam_FastTree_treeDist/allPfamsAveraged_treeDist.txt'

python3 convert2phylip10Padding.py $data_dir'/pfam_FastTree_treeDist/allPfamsAveraged_treeDist.phylip' $data_dir

if [ -f 'outfile' ]; then
    rm -r outfile
    rm -r outtree
    rm $tree_dir
    mkdir $tree_dir
fi

neighbor < neighbor_cmds.txt > screenout &

mv outfile $tree_dir'allPfamsAveraged_treeDist.outfile'
mv outtree $tree_dir'allPfamsAveraged_treeDist.outtree'

python3 mapBackTreeLeafNames.py $data_dir'allPfamsAveraged_treeDist_padded_number2bin_dic.json' $tree_dir'allPfamsAveraged_treeDist.outtree'

python3 annotateTreeParents.py $tree_dir'allPfamsAveraged_treeDist.outtree' $tree_dir

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Exec time $DIFF seconds"
