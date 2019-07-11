## Human_Gut_Tree

### git hub repository for human gut tree project

* The repository will contain a full pipeline of building a bacterial tree of life starting from contigs and or binned genomes all the way to a tree in newick format.
* This github will also include helper scripts to annotae the tree (and the steps involved in annotating the tree)
* scripts needed to generate label files used in the metaproteomic application
* scripts needed to generate label files used in metagenomic aapplication

### Prerequisite Programs needed to run the pipeline:

* python3
* required python packages:
    * biopython
    * ete3
* Frag Gene Scan (FGS)
* HMMER3
* Bowtie2
* MUSCLE
* FastTree

    
Click [here](#setting-up-required-programs-and-packages) for details about installing/setting up prerequisites.

### Buidling phylogenetic Tree by marker profile distance averaging method
In this section we will describe the steps and the code used to build a phylogenetic tree given list of contigs. The algorithm starts from a set of contigs or genomes (depending on the completeness of your datasets) to a phylogenetic tree representing these contigs. The main steps involved for this pipeline is as follows:

* [Prediction protien coding genes from contigs or genomes](#prediction-protien-coding-genes-from-contigs-or-genomes)
* [Searching for marker gene profile matches within the predicted proteins](#searching-for-marker-gene-profile-matches-within-the-predicted-proteins)
* [Extracting pfam to sequene hits and best pfam to sequence hit](#extracting-pfam-to-sequene-hits-and-best-pfam-to-sequence-hit)
* [Creating pfam multi-fasta files for multiple sequence alignments](#creating-pfam-multi-fasta-files-for-multiple-sequence-alignments)
* [Performing multiple sequence alignments using MUSCLE](#performing-multiple-sequence-alignments-using-MUSCLE)
* [Phylogenetic tree construction using FastTree](#phylogenetic-tree-construction-using-fasttree)
* [Phylogenetic tree to distance matrix](#phylogenetic-tree-to-distance-matrix)
* [Combining distance matrices](#combining-distance-matrices)
* [Building a final Tree using neighborhood joining approach](#building-a-final-tree-using-neighborhood-joining-approach)

All of these steps and the code involved are explained below, to run this entire pipeline in one command we made a script that combines all the steps. This script requires 1 mandatory command line argument which is the directory where the contigs and the genome are located in fasta format. You can also satisfy the number of cores used using the -t parameters (default is 40).

To run pipleine in one command simply use the script [constructTree.sh](buildTree/treeBuildingPipeline/constructTree.sh).

example to run script:
```
sh constructTree.sh -i ../../genomes_dir/ -t 40
```

An overview of the pipeline is summarized in this figure:
![GitHub Logo](buildTree/tree_building_pipeline.png)

### Rooting the phylogenetic tree using an outgroup
In order to root out phylogenetic tree properly and to have two branches coming from the root node, we need to define an outgroup and use it to root the tree. In our bacterial tree of life example we used a set of 7 Archaea genomes whose most specific common ancestor was used as the outgroup clade. To do this we created a script called ['rootTree_usingOutgroup.py'](buildTree/treeBuildingPipeline/rootTree_usingOutgroup.py). Note that this step requires manual intervention because we need the user to spceify the outgroup genomes based on genome data used. This script requires two arguments to run:

1) input tree file 
2) directory for the output to dump the rooted tree.

This script will also prompt the used to input the genome or bin IDs for the genomes or bins to be used as the outgroup. The user should enter the IDs without the file extensions, then the script will calculate the most recent common ancestor of this list of species entered and will use that common ancestor to be the outtgroup clade.

Example to run this script:
```
pyhton3 rooTree_usingOutgroup.py ../../data/combinedTree/allPfamsAveraged_treeDist_clean_internalNodesNamed.outtree ../../data/combinedTree/
```

### Taxonomic assignment using GTDBTK

Here we are using GTDBTK taxonomies to assign our genomes/bins with taxonomies using the least common ancestor approach. To do this we issue the following command:

```
gtdbtk classify_wf --cpus 64 --genome_dir final_genomes/ --extension fasta --out_dir data/final_genomes.classify_wf/
```
where we specified the folder where all of these genomes/bins are found in fasta format, and the output directory for the generated files to be stored plus number of CPUs used. This command will generate different taxonomical assignment files for the different kingdoms. i.e. Archaea and Bacteria in our case. These files are called 'gtdbtk.ar122.summary.tsv' and 'gtdbtk.bac120.summary.tsv' respectively.

#### Bin IDs to taxonomical mappings
After running the GTDBTK program over our genomes, we will obtain two taxonomical assignments once for bacteria and the other for archaea. We made a script that will combine these two outputs and generate a dictionary mapping, that maps between the bin IDs to the different taxonomical levels assigned to that bing or genome ID, i,e, starting to phylum all the way to species level taxonomic assignments, if that information is available. This resulting dictionary will be used later to annotate and creat labels for the tree. The script is found [here](assignTaxonomies/bins2taxonomic_assignment_GTDBTK.py). This script requires three command line arguments:

1) bacterial taxonmic assignment file
2) archaea taxonomic assignment file
3) output directory to store the this dictionary mapping between bin IDs to taxonomies.

Example to run the script:

```
python3 bins2taxonomic_assignment_GTDBTK.py ../data/final_genomes.classify_wf/gtdbtk.bac120.summary.tsv ../data/final_genomes.classify_wf/gtdbtk.ar122.summary.tsv ../data/
```

#### Bin IDs to taxonomical mappings

After obtaining bin ID to taxonomic assignments, now we can start with annotating and coloring the tree at different taxonomicla levels. To do that we made a script called [extractSpecificLevelAnnotation.py](assignTaxonomies/extractSpecificLevelAnnotation.py). This script requires three command line arguments to run:

1) First command line argument is a dictionary mapping between bin IDs to different levels of taxonomic assignments, i.e. species, genus, order etc...
2) Second command line argument is the level at which the annotation is to be extractesd i.e. phylum class etc...
3) Third input is the output diretory to store the generated files.

Example to run the script:

```
python3 extractSpecificLevelAnnotation.py ../data/allBin2taxon_dic.json phylum ../data/
```

#### Annotating internal parent nodes (LCA approach)

After having taxonomica information for the leaf nodes (i.e. out species) and after having constucted a phylogenetic tree composed of all these species, we can propagate the taxonomic information to the internal parent nodes using the tree topology and the information present at the leaf nodes by least common ancestors approach. To do this we wrote a script called ['nodes2LCA_maps.py'](assignTaxonomies/nodes2LCA_maps.py). THis script takes three arguments to runL

1) mapping dictionary file between bin IDs to taxonomies.
2) file to the final rooted phylogenetic tree containing all species.
3) directory to the output to store the node to LCA maps.

Example to run the script:

```
python3 nodes2LCA_maps.py ../data/allBin2taxon_dic.json ../data/combinedTree/allPfamsAveraged_treeDist_clean_internalNodesNamed_rooted.outtree ../data/
```


### ---------------------Detail explanation of the steps mentioned above-----------------------------

#### Prediction protien coding genes from contigs or genomes
Before we start with the gene prediction we made a small helper script that will unify the genome extensions in case your genomes come from multiple sources with different extensions, this will help to be able to split genome IDs for later. The script used for that can be found [here](buildTree/treeBuildingPipeline/unifyGenomeExtensions.py). To run this script you need to specify one command line argument which is the directory where you stored the genomes. It will overwrite the same directory by using '.fasta' as the file extension.

Example to run this script:
```
python3 unifyGenomeExtensions.py ../../sample_genomes/
```

It is essential that we have a list of genes for the genomes of interest to extract marker gene information later to build the tree. Here we use FragGeneScan (FGS) in order to predict the protein coding genes.
To predict genes we made a [script](buildTree/treeBuildingPipeline/runFGS_parallel.sh), 'runFGS_parallel.sh' that calls FragGeneScan over all the contigs specified in a particular folder. This will execute FragGeneScan in parallel treating each genome/bin independantly, it will scale as much as the specified number of cores. This script takes four command line arguments:

1) input directory containing all contigs
2) number of threads used by the program.
3) output directory for the FGS to dumpt it's predicted genes and other produced files
4) extension for the files in a folder, if your contigs have more than one extension then simply do not use this parameter

Example to run this script:
```
sh runFGS_parallel.sh -i ../../sample_genomes -t 40 -o ../../data/sample_genomes.FGS

or if extension is to be specified

sh runFGS_parallel.sh -i ../../sample_genomes -t 40 -o ../../data/sample_genomes.FGS -e .fa
```

#### Searching for marker gene profile matches within the predicted proteins
The next step after gene prediction would be to scan marker gene profiles against proteins of these genes. Since we use FGS for the task of gene prediction, it also produces the translated protein sequences of these genes. 
These produced protein sequences are used to scan (search) these marker gene profiles against them. I have predefined marker gene profiles and a precompiled hmmer database which I am including with this github. The list of pfam profiles can be found [here](buildTree/treeBuildingData/ribosomal_GTP_EFTU_pfamIDs_list.txt) and the hmmer3 database is found [here](buildTree/treeBuildingData/ribosmal_GTP_EFTU_pfam_db/)
To search for significant hits betwen sequences and these predefined profiles hmmscan function from Hmmer3 is used. To do this we made the following script [script](buildTree/treeBuildingPipeline/runHMMSCAN_parallel.sh), 'runHMMSCAN_parallel.sh'. This will call hmmscan function of hmmer, and will produce two files per scanned proteome file in this case, one is a tabular file and another is a human readable file at an output folder specified by the output parameter. This script takes five command line arguments:

1) input directory containing protein sequences in fasta format
2) number of threads used by the program.
3) output directory to dump the hmmscan output
4) directory for the precompiled hmmer database
4) extension for the files in a folder, if your fasta sequences have more than one extension then simply do not use this parameter

Example to run this script:
```
sh runHMMSCAN_parallel.sh -i ../../data/sample_genomes.FGS -t 40 -o ../../data/sample_genomes.FGS_hmmscan_out -m ../treeBuildingData/ribosmal_GTP_EFTU_pfam_db/ribosomal_GTP_EFTU_profiles.hmm

or if extension is to be specified

sh runHMMSCAN_parallel.sh -i ../../data/sample_genomes.FGS -t 40 -o ../../data/sample_genomes.FGS_hmmscan_out -m ../treeBuildingData/ribosmal_GTP_EFTU_pfam_db/ribosomal_GTP_EFTU_profiles.hmm -e .faa
```

#### Extracting pfam to sequene hits, and best pfam to sequence hit
Now that we have all the hmmscan outputs, we create a [script](buildTree/treeBuildingPipeline/extractPfamSeqHits.py) called 'extractPfamSeqHits.py', that will go over the text output of hmmscan and will parse them and extract two files from each hmmscan out file. The first one will be all the hits between the pfams and the sequenes scanned against that are significant and the second one will be for each pfam the best sequence hit between that pfam and that particualr genome. The results for this scrip are to be stored in a directory specified at command line. The script will further create two directories within the specified ourput directory. This script takes 3 command line arguments:

1) input directory for the hmmscan output files
2) number of threads used by the program
3) output directory to dump the extracted sequence hits:

example to run this script:
```
python3 extractPfamSeqHits.py ../../data/sample_genomes.FGS_hmmscan_out/ 40 ../../data/
```
#### Creating pfam multi-fasta files for multiple sequence alignments
Here we use the information from the previous script and create a multi-fasta file for each pfam and their best sequence hits with each genome. Note that some pfams might not be present in some of the genomes hence they will not have representative sequecnes from those genomes. This is done by two scripts ['binID2BestPfamSeqs.py'](buildTree/treeBuildingPipeline/binID2BestPfamSeqs.py) and ['extract_profile_sequences.py'](buildTree/treeBuildingPipeline/extract_profile_sequences.py)

The script 'binID2BestPfamSeqs.py' takes two command line arguments:
1) directory for the best pfam hits directory
2) output directory to store the files

The script 'extract_profile_sequences.py' takes two command line arguments:
1) the pfams to binns to sequences dictionary of dictionaries file
2) output directory to store the files

```
python3 binID2BestPfamSeqs.py ../../data/gene2bestpfam_hits/ ../../data/

and 

python3 extract_profile_sequences.py ../../data/pfam2bins2bestPfamSeqs_dic_of_dics.json ../../data/bin2bestPfam_seqs/
```


#### Performing multiple sequence alignments using MUSCLE
Now that best sequence hits between pfams and genomes are extracted in mulit-fasta format, MUSCLE is used to calculate one multiple sequence alignment per pfam. To do this we use two scripts, a python script calling MUSCLE which is ['get_pfam_MSA.py'](buildTree/treeBuildingPipeline/get_pfam_MSA.py) and another bash script that parallelizes this process found [here](buildTree/treeBuildingPipeline/get_MSA_parallel.sh). The bash script takes three arguments:

1) input directory for the folder containing multi-fasta files
2) number of threads to be used
3) output directory for the script to dump the multiple sequence alignments (in fasta format)

example to run the script:
```
sh get_MSA_parallel.sh -i ../../data/bin2bestPfam_seqs -t 40 -o ../../data/pfam_MSA
````

#### Phylogenetic tree construction using FastTree
After obtaining the multiple sequence alignments for each pfam, we use FastTree to construct one phylogenetic tree per pfam. To do this we use two scripts, a python script calling MUSCLE which is ['getFastTreeFromMSA.py'](buildTree/treeBuildingPipeline/getFastTreeFromMSA.py) and another bash script that parallelizes this process found [here](buildTree/treeBuildingPipeline/get_FastTree_parallel.sh). The bash script takes three arguments:

1) input directory for the folder containing multiple sequence alignments
2) number of threads to be used
3) output directory for the script to dump the phylogenetic trees (in newick format)

example to run the script:
```
sh get_FastTree_parallel.sh -i ../../data/pfam_MSA -t 20 -o ../../data/pfam_FastTree
````

#### Phylogenetic tree to distance matrix
After having created individual trees for each pfam, now we need to get pairwise species evolutionary distances from these trees. We do that by using and R package which calculates such distances based on tree branches. We use an Rscript to do that which can be found [here](buildTree/treeBuildingPipeline/cophenetic_phylo.R). This script takes two command line arguments to run:

1) the input directory for the phylogenetic trees
2) the output directory to store the calculated distance matrices

example to run the script:
```
Rscript cophenetic_phylo.R ../../data/pfam_FastTree/ ../../data/pfam_FastTree_treeDist/
```

#### Combining distance matrices
Now that we have individual distance matrices for each of the pfam profiles, it is time to do the averaging step and complie one big matrix composed of all the species, whose values are going to be the average values of all the values within the individual matrices for a particular species pair. Note that if a distance between a species pair does not exist in one of the pfam matrices, then that matrix will not be included in the averaging for that species pair. To perform this averaging we made the script called ['combTreeDistMats.py'](buildTree/treeBuildingPipeline/combTreeDistMats.py). This script takes one command line argument to run:

1) input directoru for the calculated pfam distance matrices

example to run the script:
```
python3 combTreeDistMats.py ../../data/pfam_FastTree_treeDist/
```

#### Building a final Tree using neighborhood joining approach
After combining all the individual pfam matrices into one matrix containing all against all pairwise evolutionary distances between the pariticipating species, we now have all the necesarry information to build a final tree covering all the species. To build the final tree from a distance matrix we are using a neighborhood joining approach provided by PHYLIP software. This requires a few steps to get to a tree from a distance matrix. The first step is to convert the distance matrix to phylip format, this is done by the script ['dfMat2phylip.py'](buildTree/treeBuildingPipeline/dfMat2phylip.py). This script requires one command line argumment which is the averaged distance matrix file tor run:

1) directory to the averaged pfam distance matrix file

We also need to change the leaf names to padded format since neighborhood joining program in PHYLIP requires a padding of 10 for the leaf IDs within the distance matrix. This is achieved through the script ['convert2phylip10Padding.py'](buildTree/treeBuildingPipeline/convert2phylip10Padding.py). This script requires two command line arguments:

1) path to the input matrix file created by the first script
2) output  directory to store the mapping dictionary between these arbitrary padded leaf IDs and the actual leaf IDs.

This script will also produce a file that is the command line arguments to run PHYLIP's neighbor program.

Next it is time to call PHYLIP's neighbor program through the command line to actually create a phylogenetically tree in newick format.
The final step is to map back these arbitrary IDs to the actual leaf IDs using the created dictionary. This is done through the script ['mapBackTreeLeafNames.py'](buildTree/treeBuildingPipeline/mapBackTreeLeafNames.py), which requires two command line arguments to run:

1) path to the ID mapping dictionary file
2) path to the output tree created by PHYLIP's neighbor program

A final step after creating the tree is to provide the internal parent nodes with names so that later they can be referred. To do this we wrote a script called ['annotateTreeParents.py'](buildTree/treeBuildingPipeline/annotateTreeParents.py), which requries two command line arguments to run:

1) path to the tree file
2) output directory to store the named tree.

example to run these scripts:
```
python3 dfMat2phylip.py ../../data/pfam_FastTree_treeDist/allPfamsAveraged_treeDist.txt

python3 convert2phylip10Padding.py ../../data/pfam_FastTree_treeDist/allPfamsAveraged_treeDist.phylip ../../data/

neighbor < neighbor_cmds.txt > screenout &

mv outfile ../../data/combinedTree/allPfamsAveraged_treeDist.outfile
mv outtree ../../data/combinedTree/allPfamsAveraged_treeDist.outtree

python3 mapBackTreeLeafNames.py ../../data/allPfamsAveraged_treeDist_padded_number2bin_dic.json ../../data/combinedTree/allPfamsAveraged_treeDist.outtree

python3 annotateTreeParents.py ../../data/combinedTree/allPfamsAveraged_treeDist.outtree ../../data/combinedTree/
```



## Setting up required programs and packages

#### Python3 and Anaconda3
Installing python3 through anaconda suite, will install alot ofthe python packages used and will make it easy to install additional packages through conda. Anaconda3 package could be obtained through:
```
wget https://repo.anaconda.com/archive/Anaconda3-2018.12-Linux-x86_64.sh
sh Anaconda3-2018.12-Linux-x86_64.sh
and then follow the installation instructions.
```

##### Installing biopython and ete3 python packages
```
conda install -c anaconda biopytho
conda install -c etetoolkit ete3 
```


#### Frag Gene Scan
Frag Gene Scan is used to predicted protein coding genes from the contigs. It could be downloaded from sourceforge:
```
wget https://sourceforge.net/projects/fraggenescan/files/latest/download/FragGeneScan1.31.tar.gz
tar -zxvf FragGeneScan1.31.tar.gz
add the path to FragGeneScan to $PATH variable in ~/.bashrc file
export PATH="/home/mstambou/tree_match_pipeline/tree_match_prereqs/FragGeneScan1.31:$PATH"
```

#### hmmer3
hmmer3 is used for creating marker gene database also searching for marker genes based on a predefined database, using the hmmscan functionality. hmmer could be downloaded from:
```
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar -zxvf hmmer.tar.gz
cd hmmer-3.2.1/
./configure --prefix
./configure --prefix /home/mstambou/tree_match_pipeline/tree_match_prereqs/hmmer-3.2.1
 make
 make check
 make install
 add the path to hmmer to the $PATH variable in ~/.bashrc file
 export PATH="/path/to/hmmer/hmmer-3.2.1/bin:$PATH
```
#### BOWTIE2
BOWTIE2 is used to map the reads back to the contigs, and create SAM files. BOWTIE2 could be obtained from:
```
wget https://sourceforge.net/projects/bowtie-bio/files/latest/download/bowtie2-2.3.4.3-source.zip
unzip bowtie2-2.3.4.3-source.zip
add assign the path to bowtie2 to the variable $BT2_home, in ~/.bashrc file.
export BT2_HOME="/path/to/bowtie2/bowtie2-2.3.4.3"
make
```
If make gives errors you might have to switch back to the default gcc compiler, rather than homebrew's one, and then make.

#### MUSCLE
MUSCLE is used to peform multiple sequence allignment between the sequences belonging to the same maker gene classes. MUSCLE could be obtained from:
```
wget https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
tar -zxvf muscle3.8.31_i86linux64.tar.gz
make the binary executable
chmod +x /path/to/muscle/muscle3.8.31/muscle3.8.31_i86linux64
add the path to MUSCLE to the $PATH variable in ~/.bashrc file
export PATH="/path/to/cd-hit/muscle3.8.31:$PATH"
```

#### FastTree
FastTree is used to build a phylogenetic tree from the multiple sequence allignments obtained by MUSCLE. FastTree could be obtained from:
```
wget http://www.microbesonline.org/fasttree/FastTree
chmod +x /path/to/FastTree/FastTree
add the path to FastTree to the $PATH variable in ~/.bashrc file
export PATH="/path/to/FastTree:$PATH"
```

### Prerequisite Programs needed to run the GTDBTK to annotate species:

* python2
* required python packages:
    * dendropy
    * future
    * Matplotlib
    * Numpy
    * Scipy
* Prodigal
* pplacer
* FastANI
* FastTree
* GTDBTK

Click [here](#setting-up-GTDBTK) for details about installing/setting up prerequisites.


## Setting up GTDBTK
#### Python2 and Anaconda2
Installing python3 through anaconda suite, will install alot ofthe python packages used and will make it easy to install additional packages through conda. Anaconda3 package could be obtained through:
```
wget https://repo.anaconda.com/archive/Anaconda2-2019.03-Linux-x86_64.sh
sh Anaconda2-2019.03-Linux-x86_64.sh
and then follow the installation instructions.
```
to install pyhton2 packages:
```
conda install -c bioconda dendropy
pip install future
conda install -c conda-forge matplotlib
conda install -c anaconda numpy 
conda install -c anaconda scipy 

```

#### Prodigal

to install Prodigal follow these steps:
```
git clone https://github.com/hyattpd/Prodigal.git.
cd Prodigal
make install
make install
```

#### installing pplacer

to install pplacer we need to install a number of other packages first:

```
wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip
unzip the folder
make sure to add the directory into the path variable
```

#### installing FastANI

instructions:
```
git clone https://github.com/ParBLiSS/FastANI.git
cd FastANI
./configure
make  (note that when making the binaries, this might throw an error depending on the version of your g++, make sure you have a gcc compiler version 5 or newer)
```

#### FastTree

download FastTree exectuable and add it's directory to the path variable.

#### GTDBTK

download and install GTDBTK.

data needed to be downloaded:
```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
tar xvzf gtdbtk_r89_data.tar.gz
wget https://data.ace.uq.edu.au/public/gtdbtk

```
installing the package:

pip install gtdbtk

download and unpack taxonomical information for gtdbtk to function:

wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
tar xvzf gtdbtk_r89_data.tar.gz
wget https://data.ace.uq.edu.au/public/gtdbtk

set the GTDBTK_DATA_PATH variable to point where these datasets were downloaded.
```


