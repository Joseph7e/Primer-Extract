# Primer-Extract

## Starting Database: https://www.arb-silva.de/download/archive/qiime/
SILVA_18S_99_OTUS.fasta

SILVA_18S_majority_taxonomy.tsv

## 18S Primers: 5' -> 3'
### V1-V2 primers:
 Forward: GCTTGTCTCAAAGATTAAGCC
 
 Reverse:CCTGCTGCCTTCCTTRGA
### V9 primers:
 Forward:GTACACACCGCCCGTC
 
 Reverse:TGATCCTTCTGCAGGTTCACCTAC


## Narrow reference to Metazoa
python3 extract_metazoa.py

## Extract amplicon region (run twice, once for each primer set)

### Method 1: primer-blast (be sure to rev comp reverse primers)
blastn -task blastn-short -query $concatenated_primers -subject ../metazoa.fasta -outfmt '6 qseqid qlen sseqid slen length pident qstart qend sstart send evalue bitscore' -out metazoa_v9.blast -word_size 7 -evalue 1000

### Method 2: Cutadapt
cutadapt --discard-untrimmed -g $forward_primer 99_otus_18S.fasta 2> /dev/null | cutadapt --discard-untrimmed -a $reverse_primer - 2> /dev/null > extracted_region.fasta


## Method 3: CLC in silco PCR extraction
CLC has a tool that uses BLAST, similiar to above. Results are consitent across the three methods.


## Parse BLAST (if necessary) and divide extracted sequence into a per phylum fasta file.
python3 primerRegionExtractor.py

## Align and construct distance matrix with clustalo, once per phylum
clustalo -i $extracted_region_per_phylum -o Annelida_regions_V9.align --distmat-out=Annelida_distance_V9.txt --full
### Run on all fastas
./auto_run_alignment.sh <fasta1> <fasta2> ...

### Parse all distance matrixes and construct output figures (Original uses jupyter notebook)
#### On one file:
python3 calculate_average_distance.py
#### on all files:
python3 Distance_Matrix_phylum_based.py

### Convert jupyter notebook to regular python script
jupyter nbconvert --to script Distance_Matrix_phylum_based.ipynb
