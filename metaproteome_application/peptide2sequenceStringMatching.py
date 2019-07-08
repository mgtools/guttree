"""
script that will perform an exact string matching between the peptides
and protein sequences from the genomes and will rerturn a table
reporting these matches.

this script requires three command line arguments
the first argument being the file that has the peptides to be searched in fasta format
these peptides could be obtained by any peptide searching methods such as MSGF+ or some other way
second command line argument is the file containing all the proteins coming fromall the genomes
I pre-compiled a file by predicting from all the contigs the protein coding genes and then translating these
genes into their relative amino acid sequences. If you have not done so you need to precomple such a file
before runing this script
the third command line argument is the directory for the output file:
"""

from Bio import SeqIO
import sys

#peptides_seq_f = 'unipeptSample7Peptides2myGenomes/mpa_results_peptide_list.fasta'
#all_prot_seqs_f = 'allBinsProtSeqs.faa'

if len(sys.argv) != 4:
    print('please enter 3 command line arguments to run this script, example to run script i.e. python3 peptides2sequenceStringMatching.py dir/2/peptides/file dir/to/allProteins/file/ dir/to/output/directory')

else:
    peptides_seq_f = sys.argv[1]
    all_prot_seqs_f = sys.argv[2]
    out_dir = sys.argv[3]
    out_f_name = peptides_seq_f.rsplit('/', 1)[-1]
    
    peptide_seqs = SeqIO.to_dict(SeqIO.parse(peptides_seq_f, 'fasta'))
    print('reading all proteins into memory')
    all_prot_seqs = list(SeqIO.parse(all_prot_seqs_f, 'fasta'))
    
    peptide_seq_list = [(peptide_id ,str(peptide_seqs[peptide_id].seq)) for peptide_id in peptide_seqs]
    
    new2oldPeptide_dic = dict()
    new_peptide_list = list()
    
    print('converting peptides I 2 L')
    for peptide in peptide_seq_list:
        new_peptide = peptide[1].replace('I', 'L')
        new_tuple = (peptide[0], new_peptide)
        new2oldPeptide_dic[new_peptide] = peptide[1]
        new_peptide_list.append(new_tuple)
    
    print('converting protein seqs I 2 L')
    for seq in all_prot_seqs:
        seq.seq = str(seq.seq).replace('I', 'L')
    
    print('searching peptides against the proteins')
    with open(out_dir+out_f_name, 'w') as out_f:
        out_f.write('peptide_ID\tpeptide\tprotein_match_ID\n')
        for i, peptide in enumerate(new_peptide_list):
            print(i, peptide[0], peptide[1])
            for seq in all_prot_seqs:
                prot_seq = str(seq.seq)
                prot_id = str(seq.id)
                if peptide[1] in prot_seq:
                    out_f.write(str(peptide[0])+'\t'+str(new2oldPeptide_dic[peptide[1]])+'\t'+str(prot_id)+'\n')