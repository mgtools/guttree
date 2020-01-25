#!/bin/bash
#
#script that will be used to bootstrap the multiple sequence alignment used to construct a tree
#and calculate the branch confidence intervals.
#

START=$(date +%s)

data_dir='../data/'
bootstrap_dir='../data/bootstrap/'

while getopts i:t: option
do
case "${option}"
in
    i) tree_f=${OPTARG};;
    t) n_cores=${OPTARG};;
esac
done


if [ ! -d $bootstrap_dir ]; then
    mkdir -p $bootstrap_dir
fi

cp -r $data_dir'supermat_clean_columns_clean_rows_FastTree.tree' $bootstrap_dir
cp -r $data_dir'supermat_clean_columns_clean_rows.fasta' $bootstrap_dir

python3 seqID2num.py $bootstrap_dir'supermat_clean_columns_clean_rows.fasta'

python3 fasta2phylip.py $bootstrap_dir'supermat_clean_columns_clean_rows_paddedIDs.fasta'

rm -r outfile

seqboot < seqboot_cmds.txt

mv outfile supermat_clean_columns_clean_rows_paddedIDs_100bootstraps.phylip

mv supermat_clean_columns_clean_rows_paddedIDs_100bootstraps.phylip $bootstrap_dir

python3 splitAlignments.py $bootstrap_dir"supermat_clean_columns_clean_rows_paddedIDs_100bootstraps.phylip" $bootstrap_dir"alignment_out/"

sh get_FastTree_parallel.sh -i $bootstrap_dir"alignment_out/" -t $n_cores -o $bootstrap_dir"FastTree_out/"

cat $bootstrap_dir"FastTree_out/"*.tree > $bootstrap_dir"supermat_clean_columns_clean_rows_paddedIDs_100bootstraps_FastTree.tree"

python3 mapBackLeafNames.py $bootstrap_dir"supermat_clean_columns_clean_rows_paddedIDs_100bootstraps_FastTree.tree" $bootstrap_dir"supermat_clean_columns_clean_rows_paddedNum2ID_dic.json"

perl CompareToBootstrap.pl -tree $bootstrap_dir"supermat_clean_columns_clean_rows_FastTree.tree" -boot $bootstrap_dir"supermat_clean_columns_clean_rows_originalIDs_100bootstraps_FastTree.tree" > $bootstrap_dir"supermat_clean_columns_clean_rows_FastTree_withSupports.tree"
