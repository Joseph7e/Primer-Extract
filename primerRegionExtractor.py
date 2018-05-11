import sys
import os

# The purpose of this program is to cut out a primer region from a fasta (if necessary) and output a per phylum FASTA file of regions

# V9 forward primer - GTACACACCGCCCGTC
#cutadapt --discard-untrimmed -g GTACACACCGCCCGTC ../metazoa.fasta > trimmed_metazoa_V9.fasta
# V1-V2 reverse primer - CCTGCTGCCTTCCTTRGA
# revseq primers_V1_V2.fasta --> TCYAAGGAAGGCAGCAGG
# cutadapt --discard-untrimmed -a TCYAAGGAAGGCAGCAGG ../metazoa.fasta > trimmed_metazoa_V2.fasta



from Bio import SeqIO
# input a fasta file, a forward read, and a reverse read
# output fasta_file with only amplified region and a log of each sequence specificity


working_dir = "/mnt/lustre/hcgs/joseph7e/program_PrimerRegionExtractor/extract_from_blast/"

taxonomy_file = working_dir + "../majority_taxonomy_all_levels.txt"
#interested_taxa = sys.argv[2]
# blast_results = working_dir + "V9.blast"

def OutputMatching_regions(dict1, dict2):
    output1_handle = open('V9_regions_' + interested_taxa +  '.fasta','w')
    output2_handle = open('V2_regions_' + interested_taxa +  '.fasta', 'w')
    for k, v in dict1.items():
        #if tax_dict[k] == interested_taxa:
        if k  in dict2.keys():
            output1_handle.writelines('\n>' + k + '\n' + v)
            output2_handle.writelines('\n>' + k + '\n' + dict2[k])


# populate taxonomy dictionary
tax_dict = {}
for line in open(taxonomy_file):
    elements = line.rstrip().split('\t')
    seq = elements[0]
    phylum = elements[1].split(';')[6]
    if phylum != 'Ambiguous_taxa':
        phylum = phylum.split('__')[1]
    phylum = phylum.replace(' ', '_')
    tax_dict[seq] = phylum

# Populate ssequence dictionary
fasta_file_V2 = "trimmed_metazoa_V2.fasta" # fasta file
fasta_file_V9 = "trimmed_metazoa_V9.fasta" # fasta file
seq_dict_V9 = SeqIO.to_dict(SeqIO.parse(fasta_file_V9, "fasta"))
seq_dict_V2 = SeqIO.to_dict(SeqIO.parse(fasta_file_V2, "fasta"))

# ensure extracted region makes length requirement and that each output will have same ref sequences
matching_list = []
V9_length = 99
V2_length = 299
for seqid in seq_dict_V2.keys():
    if seqid in seq_dict_V9.keys():
        if len(seq_dict_V9[seqid]) > V9_length and len(seq_dict_V2[seqid]) > V2_length:
            matching_list.append(seqid)


# Construct output
outdir_V2 = "V2_regions/"
outdir_V9 = "V9_regions/"
os.system("rm V2_regions/* V9_regions/*")

handles_V2 = []
handles_V9 = []

for seqid in matching_list: #
    # deal with V2
    sequence = seq_dict_V2[seqid]
    phylum = tax_dict[seqid]
    out_name = outdir_V2 + phylum + '.fasta'
    if out_name not in handles_V2:
        handles_V2.append(out_name)
    output_handle = open(out_name, 'a')
    output_handle.writelines('\n>' + seqid + '\n' + sequence)

    # deal with V9
    sequence = seq_dict_V9[seqid]
    phylum = tax_dict[seqid]
    out_name = outdir_V9 + phylum + '.fasta'
    if out_name not in handles_V9:
        handles_V9.append(out_name)
    output_handle = open(out_name, 'a')
    output_handle.writelines('\n>' + seqid + '\n' + sequence)


# Definitions if you use BLAST or CLC methods
def getCutSitesfromBLAST(blast_results, direction, keep_side):
    cut_dict = {} # seq_id: [cut_site, cut_side] ex. HP642068.311: [320, left]
    # keep side right for example will remove everything to the left of site
    for line in open(blast_results):
        qseqid, qlen, sseqid, slen, length, pident, qstart, qend, sstart, send, evalue, bitscore = line.rstrip().split()
        if keep_side ==  "right":
            cut_dict[sseqid] = [sstart,"right"]
        else:
            cut_dict[sseqid] = [send, "left"]
    return cut_dict


def getCutSitesfromCLC(clc_file, direction, cut_side):
    out_dict = {}
    for line in open(clc_file):
        elements = line.rstrip().split('\t')
        if elements[1] != "Name":
            seq = elements[0]
            orientation = elements[3]
            region = elements[4]
            if orientation == direction:
                if "complement" in region:
                    start = region.split('..')[0].split('(')[1]
                    end = region.split('..')[1].replace(')','')
                else:
                    start, end = region.split('..')

                if cut_side == "right":
                    out_dict[seq] = [int(start),cut_side]
                else:
                    out_dict[seq] = [int(end), cut_side]
    return out_dict

def cutSequence(fasta_dict, cutsites_dict):
    new_seqs = {}
    for seqid, data in cutsites_dict.items():
        original_sequence = fasta_dict[seqid]
        cut_site = data[0]
        side = data[1]
        if side == "right": # trim all data on the left of cut site
            new_sequence = original_sequence[cut_site:]
            new_seqs[seqid] = new_sequence
        else: # trim all the data on the right side of the cut site
            new_sequence = original_sequence[:cut_site]
            new_seqs[seqid] = new_sequence
    return (new_seqs)


def OutputMatching_regions(dict1, dict2):
    output1_handle = open('V9_regions_' + interested_taxa +  '.fasta','w')
    output2_handle = open('V2_regions_' + interested_taxa +  '.fasta', 'w')
    for k, v in dict1.items():
        #if tax_dict[k] == interested_taxa:
        if k  in dict2.keys():
            output1_handle.writelines('\n>' + k + '\n' + v)
            output2_handle.writelines('\n>' + k + '\n' + dict2[k])

# Commented out code, uncomment to use other methods, requires tweaking input

# v9_cut_sites_dict = getCutSites(V9_file, "fwd", "right")
# v2_cut_sites_dict = getCutSites(V2_file, "rev", "left")
# #v2S_cut_sites_dict = getCutSites(V2_sanger_file, "rev", "left")
#
# v9_seqs = cutSequence(sequence_dict, v9_cut_sites_dict)
# v2_seqs = cutSequence(sequence_dict, v2_cut_sites_dict)
# #v2S_seqs = cutSequence(sequence_dict, v2S_cut_sites_dict)
#


# OutputMatching_regions(v9_seqs, v2_seqs)


#blastn -query primers.fasta -db ~/databases/qiime_18S/99_otus_18S.fasta -outfmt '6 qseqid qlen sseqid slen length pident qstart qend sstart send evalue bitscore' -out three -word_size 7 -evalue 1000