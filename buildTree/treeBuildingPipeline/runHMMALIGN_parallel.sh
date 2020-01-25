#!/bin/bash 


while getopts i:t:o:e: option
do
case "${option}"
in
    i) in_dir=${OPTARG};;
    t) n_threads=${OPTARG};;
    o) out_dir=${OPTARG};;
    e) extension=${OPTARG};;
esac
done


if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi



files=( $in_dir*$extension)
n="${#files[@]}"
n="$(($n-1))"
#printf '%s\n' "${files[@]}"
echo $n
for y in $(eval echo "{0..$n}");
do
   printf "%s\0"  "${files[$y]}"; done| xargs -0 -I @@ -P $n_threads sh runHMMALIGNcommand.sh -i @@ -o $out_dir

