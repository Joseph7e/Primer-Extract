
# coding: utf-8

# In[4]:

# output from: calculate_average_distance.py
import os, statistics

working_dir = "/mnt/lustre/hcgs/joseph7e/program_PrimerRegionExtractor/extract_from_blast/"

distance_path = working_dir + "distance_matrixes/"
matrix_files = [distance_path + x for x in os.listdir(distance_path)]

target_phyla= ["Annelida", "Gastrotricha", "Mollusca", "Nematoda","Nemertea","Platyhelminthes","Xenacoelomorpha", "Tardigrada"]
all_phyla = []
primer_list = ["V9", "V1-V2"]
# parse distance matrixes
phylum_lookup = {} # phylum: [int1,int2]

for file in matrix_files:
    phylum = file.split('.fasta')[0].split('/')[-1]
    all_phyla.append(phylum)
    if phylum in target_phyla:
        seq_count = 0
        distances = []
        primer = 'V1-V2'
        if file.endswith('.V9'):
            primer = 'V9'
        for line in open(file):
            elements = line.rstrip().split(' ')
            if len(elements) != 1:
                seq_count += 1
                flag = True
                for num in elements[1:]:
                    try:
                        num = float(num)
                        if num > 0 and flag:
                            num = 1 - num
                            distances.append(num)
                        else:
                            flag = False
                    except:
                        #print ("Offending Value:", num)
                        why = True
        phylum_lookup[phylum+'   |   '+primer] = [seq_count, distances]
print (sorted(all_phyla))
print ("-----------------------------------------------")


# In[22]:




# In[5]:

headers = ["Phylum","Primer", "Number_of_Seqs","Total_comparisons","Average_Distance","STDV"]
print ('\t'.join(headers))

data_labels = []
data_matrix = []
data_labels.append('ALL_V9')
data_labels.append('ALL_V1-V2')
data_matrix.append([])
data_matrix.append([])
for phylum in target_phyla:
    two_primer_averages = []
    for primer in primer_list:
        key = phylum + '   |   ' + primer
        data_labels.append(key)
        distances = phylum_lookup[key][1]
        seq_count = phylum_lookup[key][0]
        data_matrix.append(distances)
        if primer == "V9":
            data_matrix[0].extend(distances)
        else:
             data_matrix[1].extend(distances)
        data = [phylum, primer,seq_count, len(distances), statistics.mean(distances), statistics.stdev(distances)]
        two_primer_averages.append(statistics.mean(distances))
        print ('\t'.join([str(x) for x in data]))
    print ('Primer Differences: ', two_primer_averages[1] - two_primer_averages[0], '\n' )


# In[3]:

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
def parseData(data_table):
    
    for line in open(data_table):
        elements = line.rstrip().split(',')
        phylum = elements[0].split('.fasta')[0]
        

def boxPlot(x_data_labels,all_data,title,xlabel,ylabel, v_amount):
    # all data =  list of lists
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    bp = ax.boxplot(all_data)
    fig.set_figheight(11)
    fig.set_figwidth(30)
    ax.set_xticklabels(x_data_labels,rotation=90, fontsize=20)
    ax.set_xticks(np.arange(len(x_data_labels)), minor=True)
    #add title and labels

    vertical_points = np.arange(len(x_data_labels))
    for v in vertical_points:
        if v%v_amount == 0:
            plt.axvline(x=v+.5) # adds vertical lines
    fig.suptitle(title, fontsize=40)
    #plt.xlabel(xlabel, fontsize=30)
    plt.ylabel(ylabel, fontsize=30)
    #plt.gcf().subplots_adjust(bottom=0.35)
    #plt.ylim((-100,100))
    #plt.savefig('Genes.png',bbox_inches='tight')
    plt.show() 

print ("_______________________________")

boxPlot(data_labels, data_matrix, 'Primer Region Comparison: "V9"  | "V1-V2"', "phylum", "Percent Identity", 2)

