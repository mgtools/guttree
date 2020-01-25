#!/bin/bash 

while getopts i:t:o:e: option
do
case "${option}"
in
    i) in_dir=${OPTARG};;
    o) out_dir=${OPTARG};;
    e) extension=${OPTARG};;
esac
done

if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

for file in $(ls $in_dir*$extension)
do
        echo $file
        echo $extension
        #fname=$(echo $file | cut -d'/' -f2)                                                                                                                                     
        fname=$(echo $file | rev | cut -d'/' -f1 | rev)
        profile_name=$(echo $fname | cut -d'.' -f1)
        echo $fname
	echo $profile_name

	echo "'hmmalign' $in_dir$fname 'bin2bestPfam_seqs'$profile_name'_hits.faa' > $out_dir$profile_name'_MSA.stockholm'"
	hmmalign $in_dir$fname 'bin2bestPfam_seqs/'$profile_name'_hits.faa' > $out_dir$profile_name'_MSA.stockholm'

done

