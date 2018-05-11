

for fasta in $@
do

	clustalo -i $fasta -o $fasta.regions_V2.align --distmat-out=$fasta.distance_V2.txt --full

done
