#!/bin/bash

while getopts i:t:o: option
do
case "${option}"
in
    i) in_dir=${OPTARG};;
    t) n_threads=${OPTARG};;
    o) out_dir=${OPTARG};;
esac
done



msa_files=( $in_dir*)
echo $msa_files

n="${#msa_files[@]}"

for ((i=1; i<${#msa_files[@]}; i++)); do
echo "${msa_files[$i]}"
done

for y in $(eval echo "{0..$n}");
do 
    printf "%s\0"  "${msa_files[$y]}"; done| xargs -0 -I @@ -P $n_threads python3 get_pfam_MSA.py  @@ $out_dir
