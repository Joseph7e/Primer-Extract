# Primer-Extract

## Starting Database: https://www.arb-silva.de/download/archive/qiime/
SILVA_18S_99_OTUS.fasta
SILVA_18S_majority_taxonomy.tsv

## 18S Primers: 5' -> 3'
### V1-V2 primers:
 Forward: 
 Reverse:
### V9 primers:
 Forward:
 Reverse:

## Extract amplicon region

### Method 1: primer-blast
blastn -task blastn-short -query $concatenated_primers -subject ../metazoa.fasta -outfmt '6 qseqid qlen sseqid slen length pident qstart qend sstart send evalue bitscore' -out metazoa_v9.blast -word_size 7 -evalue 1000


### Method 2: Cutadapt
cutadapt --discard-untrimmed -g $forward_primer 99_otus_18S.fasta 2> /dev/null | cutadapt --discard-untrimmed -a $reverse_primer - 2> /dev/null > output



### Align and construct distance matrix with clustalo, once per phylum
clustalo -i $extracted_region_per_phylum -o Annelida_regions_V9.align --distmat-out=Annelida_distance_V9.txt --full
