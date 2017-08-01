
# This is the first part of the programme, it integrates mutations with their coding sequences and translates the sequence into codons.
# This section takes a fasta format file (here 'results.txt') and makes a dictionary out of it ('genedict') so an ENST label corresponds to a coding sequence.
# It then makes a dictionary that has codons correpsonding to nucleotide triplets ('codon_list').
# Next, it uses the mutation dataset to makes a list('mut_ind'), each list item contains an ENST label, the position of the mutation and what nucleotides are inserted.
# Next, it loops through 'genedict' looking for ENSTs that match an item in 'mut_ind'. 
# When it finds a match, the nucleotides are added at the position specified in that same list item, this gives a new coding sequence ('geneCode').
# 'geneCode' is then transcribed into mRNA, giving a string called 'RNA'.
# 'RNA' is then translated into codons using 'codon_list'.
# Finally, it loops through 'codon_list' to find the first stop codon, length of 'codon_list'.
# The final output is ENST, gene name(?), position of stop codon ('stopcod'), difference between mutation position and stop codon ('Diff' - measured in nucleotide number), percentage of the 
# coding sequence translated ('Perc') and the nucleotide triplet prior to the mutation position. - for every mutation.

import os
import csv
import re
import gzip
from collections import defaultdict
import numpy as np
# os.chdir("/mnt/nfs2/biochem/hb280")          
f = open('practiceresults.txt', 'r')                                  
genedict = defaultdict(str)                                        
gene_name = ''
for line in f:                                                      
    if line.startswith('>'):                 
        gene_name = line[17:32]                                  
        continue
    genedict[gene_name] += line.strip()     

                                             
codon_list = {'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'X', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
              'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
              'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
              'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
              'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
              'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
              'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S', 'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
              'UAC': 'Y', 'UAU': 'Y', 'UAA': 'Z', 'UAG': 'Z','UGC': 'C', 'UGU': 'C', 'UGA': 'Z', 'UGG': 'W'}     
Gene_mut_pos = {}
AA_mut_pos = {}


mut_ind = []                                                                                   
with open('in_frameshift.txt', 'r') as indr:                          			       
    for rowInd in csv.reader(indr, dialect = 'excel-tab'):              		       	
        geneIde = rowInd[1]                                                                    
        mutationInd = list((map(int,(re.findall("\d+", rowInd[17]. partition("_")[0])))))      
        mutIns = ''									       
        for char in rowInd[17]:                                                                
            if char.isupper() == True:
                mutIns += char
            elif char.isupper == False:
                mutIns += ''
        mut_ind += [[rowInd[4], geneIde, mutationInd, mutIns]]                                 

for genam, genCode in genedict.items():
    for a in mut_ind:
        geneCode = ''
        if (genam == a[1]) and (len(genCode) - a[2][0] > 0):                                         
            Gene_mut_pos[genam] = a[2][0]
            for p in range(len(genCode)):
                if p == a[2][0]:                                 
                    geneCode += (genCode[p] + a[3])
                elif p != a[2][0]:                                  
                    geneCode += genCode[p]
            Rna = ''
            for m in geneCode:                                          
                if m == 'T':
                    Rna += 'U'
                elif m != 'T':
                    Rna += m
            CodonInd = ''
            for n in range(0, len(Rna), 3):                             
                Codons = (Rna[n: n + 3])                                
                if Codons in codon_list:                                
                    CodonInd += codon_list[Codons]
                elif Codons not in codon_list:
                    CodonInd += ''                                      
            if 'Z' in CodonInd:
                stopCod = CodonInd.index('Z')
            if a[2][0]%3 == 0:
                prevNuc = [str(genCode[a[2][0] -4 : a[2][0] - 1])]
            elif a[2][0]%3 == 1:
                prevNuc = [str(genCode[a[2][0] -5 : a[2][0] -2]) + str(genCode[a[2][0]-1])]
            elif a[2][0]%3 == 2:
                prevNuc = [str(genCode[a[2][0] -6 : a[2][0] -3]) + str(genCode[a[2][0]-2]) + str (genCode[a[2][0]-1])]
            if len(CodonInd) > 0 and (stopCod > 0):
                Perc = (stopCod / len(CodonInd)) * 100
            Diff = (stopCod * 3) - a[2][0]
            #AA_mut_pos[a[0]] = [genam, stopCod, prevNuc,Perc, Diff]     
            print(a[0], genam, stopCod, prevNuc, Diff)
            a = []