#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  1 17:18:37 2025

@author: a549cell
"""

from Bio import SeqIO
import pandas as pd
import re
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
from Bio import pairwise2
from Bio.Align.Applications import MuscleCommandline
import subprocess
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
import scipy.stats as stats
from scipy.stats import norm






def compute_gc(seq):
    
    
    return (seq.count('g') + seq.count('c')) / len(seq) * 100
def parse_fasta_header(header):
    header=header.lstrip('>')
    pattern=r'^fig\|([^|)]+)\|\s*(.*?)\s*\[(.*?)\]$'
    match=re.match(pattern,header)
    if match:
        sample_id=match.group(1)
        product=match.group(2)
        organism_info=match.group(3)
        parts=[x.strip() for x in organism_info.split('|')]
        organism=parts[0] 
        organism_id=parts[1] 
        return sample_id,product,organism,organism_id
    else:
        return None ,header,None,None
    
def fasta_parsing(fasta_file2024):
    records=[]
    for record in SeqIO.parse(fasta_file2024,'fasta'):
        
        sample_id,product,organism,organsim_id=parse_fasta_header(record.description)
        sequence=str(record.seq)
        records.append({'sample id':sample_id,
                        'product':product,
                        'organism':organism,
                        'organism_id':organsim_id,
                        'sequence':sequence,
                        'length':len(sequence),
                        'gc_content':compute_gc(sequence)})
    return records



"""
def filter_record(record):
    m={len(str(record.seq))<500}
    if len(str(record.seq)) in m:
        return False
    
input_file='BVBRC_genome_feature1-6.fasta '
output_file='DNA_modified.fasta'
filtered_records=(record for record in SeqIO.parse(input_file,'fasta') if filter_record(record))
count=SeqIO.write(filtered_records,output_file,'fasta')
print()

"""

fasta_file_amino2024='BVBRC_amino_acids_feature-7.fasta'
data=fasta_parsing(fasta_file_amino2024)
df1=pd.DataFrame(data)
print(df1)

scores=[]
bins=[0,100,300,500,700,900,1100,1500,1800,2100,2500]


for i in range(len(df)):
    
    row=df.iloc[i]
    scores.append(row['length'])
    print(row['sample id'],row['sequence'])
    
plt.hist(scores,bins,histtype='bar')
plt.xlabel('Base pairs')
plt.ylabel('The number of each lengthes')
plt.show()
   
df_mean=df['length'].mean()
print('The mean of length:',df_mean)
df_median=df['length'].median()
print('The median of length:',df_median)
counts=Counter(scores)
print(counts)
df_standard=df['length'].std()
print('The standard deviation of length:',df_standard)
df_sample_variance=df['length'].var()
print('The variance of sample:',df_sample_variance)
df_population_parameter=df['length'].var(ddof=0)
print('The parameter of population:',df_population_parameter)






sem_data=[]



data1=df['length']
for i in data1:
    sem_data.append(i)

h=len(sem_data)
print(sem_data)


pop_std=np.std(sem_data,ddof=0)
print('The population of standard deviation:',pop_std)
sample_std=np.std(sem_data,ddof=1)
print('The population of sample deviation:',sample_std)




mean=np.mean(sem_data)
confidence_level=0.99
sem = stats.sem(sem_data)
print(sem)

ci=stats.t.interval(confidence_level,df=h-1,loc=mean,scale=sem)

margin_of_error=(ci[1]-ci[0])/2

print('The margin of error:',margin_of_error)
print(f'{confidence_level*100}% The confidence level:{ci}' )


q1=df['length'].quantile(0.25)
q3=df['length'].quantile(0.75)
iqr=q3-q1
mini=q1-1.5*iqr
maxi=q3-1.5*iqr
filtered_df=df[(df['length']>=mini)&(df['length']<=maxi)]
print(filtered_df)






s=[]
r=df['sequence']
for i in r:
    s.append(i)
        
h1n1_standard='ggaaaacaaaagcaacaaaaatgaaggcaatactagtagttctgctatatacatttgcaaccgcaaatgcagacacattatgtataggttatcatgcgaacaattcaacagacactgtagacacagtactagaaaagaatgtaacagtaacacactctgttaaccttctagaagacaagcataacgggaaactatgcaaactaagaggggtagccccattgcatttgggtaaatgtaacattgctggctggatcctgggaaatccagagtgtgaatcactctccacagcaagctcatggtcctacattgtgaaacaccagttcagacaatggaacgtgttacccaggagatttcatcgattatgaggagctaagagagcaattgagctcagtgtcatcatttgaaaggtttgagatattccccaagacaagttcatggcccaatcatgactcgaacaaaggtgtaacggcagcatgtcctcatgctggagcaaaaagcttctacaaaaatttaatatggctagttaaaaaaggaaattcatacccaaagctcagcaaatcctacattaatgataaagggaaagaagtcctcgtgctatggggcattcaccatccatctactagtgctgaccaacaaagtctctatcagaatgcagatgcatatgtttttgtggggtcatcaagatacagcaagaagttcaagccggaaatagcaataagacccaaagtgagggatcgagaagggagaatgaactattactggacactagtagagccgggagacaaaataacattcgaagcaactggaaatctagtggtaccgagatatgcattcgcaatggaaagaaatgctggatctggtattatcatttcagatacaccagtccacgattgcaatacaacttgtcaaacacccaagggtgctataaacaccagcctcccatttcagaatatacatccgatcacaattggaaaatgtccaaaatatgtaaaaagcacaaaattgagactggccacaggattgaggaatatcccgtctattcaatctagaggcctatttggggccattgccggtttcattgaaggggggtggacagggatggtagatggatggtacggttatcaccatcaaaatgagcaggggtcaggatatgcagccgacctgaagagcacacagaatgccattgacgagattactaacaaagtaaattctgttattgaaaagatgaatacacagttcacagcagtaggtaaagagttcaaccacctggaaaaaagaatagagaatttaaataaaaaagttgatgatggtttcctggacatttggacttacaatgccgaactgttggttctattggaaaatgaaagaactttggactaccacgattc aaatgtgaagaacttatatgaaaaggtaagaagccagctaaaaaacaatgccaaggaaattggaaacggctgctttgaattttaccacaaatgcgataacacgtgcatggaaagtgtcaaaaatgggacttatgactacccaaaatactcagaggaagcaaaattaaacagagaagaaatagatggggtaaagctggaatcaacaaggatttaccagattttggcgatctattcaactgtcgccagttcattggtactggtagtctccctgggggcaatcagtttctggatgtgctctaatgggtctctacagtgtagaatatgtatttaacattaggatttcagaagcatgagaaaaacac'

        

sequence1=np.random.randint(0,len(s))

print('The first two sequences will be aligned:','\n','sequence1:',s[sequence1],'\n','H1N1 standard:',h1n1_standard)
length_h1n1=Seq(h1n1_standard)
length1=Seq(s[sequence1])
alignments=pairwise2.align.globalms(length_h1n1,length1,1,-1,-10,-0.5)
best_alignment=alignments[0]


print('The best alignment1:')
print(format_alignment(*best_alignment))

a=[]
m=df['sequence']
for i in m:
    a.append(i)



        


sequence3=np.random.randint(0,len(a)+1)

print('The second two sequences will be aligned:','\n','sequence3:',sequence3,'\n','H1N1 standard:',h1n1_standard)
length_h1n1=Seq(h1n1_standard)
length2=Seq(a[sequence3])
alignments2=pairwise2.align.globalms(length_h1n1,length2,1,-1,-10,-0.5)
best_alignment2=alignments2[0]


print('The best alignment2:')
print(format_alignment(*best_alignment))

a3=[]
p=df['sequence']
for i in p:
    a3.append(i)

sequence5=np.random.randint(0,len(a3)+1)


print('The third two sequences will be aligned:','\n','sequence5:',a3[sequence5],'\n','H1N1 standard:',h1n1_standard)
length_h1n1=Seq(h1n1_standard)
length3=Seq(a3[sequence5])
alignments3=pairwise2.align.globalms(length_h1n1,length3,1,-1,-10,-0.5)
best_alignment3=alignments3[0]
print(best_alignment3)















    