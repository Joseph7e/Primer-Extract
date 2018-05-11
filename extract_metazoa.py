import sys
from Bio import SeqIO



working_dir = "/mnt/lustre/hcgs/joseph7e/program_PrimerRegionExtractor/extract_regions_from_CLC/"
fasta_file = working_dir + "../99_otus_18S.fasta"
#forward_primer = sys.argv[2]
#reverse_primer = sys.argv[3]
taxonomy_file = working_dir + "../majority_taxonomy_all_levels.txt"
# blast_results = working_dir + "V9.blast"

tax_dict = {}
for line in open(taxonomy_file):
    elements = line.rstrip().split('\t')
    seq = elements[0]
    phylum = elements[1].split(';')[5]
    if "Metazoa" in elements[1]:
        tax_dict[seq] = phylum

sequence_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

for k, v in tax_dict.items():
    print ('>' + k +'\n' + sequence_dict[k].seq)
