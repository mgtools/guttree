#!/bin/bash                                                                                                  


while getopts i:o:d:f: option
do
case "${option}"
in
    i) in_f=${OPTARG};;
    o) out_dir=${OPTARG};;
    d) in_dir=${OPTARG};;
    f) data_dir=${OPTARG};;
esac
done

#in_dir='../treeBuildingData/bac_120_markers/'
#data_dir='../../data/'
#echo $in_f
fname=$(echo $in_f | cut -d'/' -f4)                                            
profile_name=$(echo $fname | cut -d'.' -f1)
#echo $fname
#echo $profile_name
#echo "outdir"
#echo $out_dir

echo "'hmmalign' $in_dir$fname $data_dir'bin2bestPfam_seqs'$profile_name'_hits.faa' > $out_dir$profile_name'_MSA.stockholm'"

hmmalign $in_dir$fname $data_dir'bin2bestPfam_seqs/'$profile_name'_hits.faa' > $out_dir$profile_name'_MSA.stockholm'
