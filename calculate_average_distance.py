import sys
import statistics
# generate distance matrix with this command
# clustalo -i ../V2_regions/Nemertea.fasta -o Nemertea_regions.align --distmat-out=Nemertea_distance_V2.txt --full




list = sys.argv[1:]
print ("File,Number_of_Seqs,Total_comparisons,Average_Distance,STDV")
for file in list:
    seq_count = 0
    distances = []
    for line in open(file):
        elements = line.rstrip().split(' ')
        if len(elements) != 1:
            seq_count += 1
            flag = True
            for num in elements[1:]:
                try:
                    num = float(num)
                    if num > 0 and flag:
                        distances.append(num)
                    else:
                        flag = False
                except:
                    #print ("Offending Value:", num)
                    why = True
                #sys.exit()
    data = [file, seq_count, len(distances), statistics.mean(distances), statistics.stdev(distances)]
    print (','.join([str(x) for x in data]))
