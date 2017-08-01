import os
import csv
import re
import gzip
from collections import defaultdict
import numpy as np
# os.chdir("/mnt/nfs2/biochem/hb280")          #change directory to location of input files
f = open('genelistseq.txt', 'r')                                   #open the fasta file containg the coding sequences (cds)
genedict = defaultdict(str)                                         #set up a dictionary
gene_name = ''
for line in f:                                                      #set the ENST as the dictionary key and the associated cds as the value
    if line.startswith('>'):
        gene_name = line[17:32]                                     # change index for del vs. ind
        continue
    genedict[gene_name] += line.strip()
                                             #make a dictionary with mRNA codons corresponding to amino acids, X = start codon, Z = stop codon
#codon_list = {'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'X', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
#              'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
#              'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
#              'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
#              'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
#              'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
#              'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S', 'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
#              'UAC': 'Y', 'UAU': 'Y', 'UAA': 'Z', 'UAG': 'Z','UGC': 'C', 'UGU': 'C', 'UGA': 'Z', 'UGG': 'W'}

codon_list = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
              'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
              'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
              'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
              'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
              'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
              'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
              'TAC': 'Y', 'TAT': 'Y', 'TAA': 'Z', 'TAG': 'Z','TGC': 'C', 'TGT': 'C', 'TGA': 'Z', 'TGG': 'W'}

gene_mut_pos = {}
aa_mut_pos = {}
mut_dict = []
#output = open('delpostnuc.txt', 'w')
with open('def_frame', 'r') as del_r:                #open the mutation data
    for row in csv.reader(del_r, dialect='excel-tab'):          #loop through each mutation
        gene_ide = row[1]                                       #find the ENST in the mutation data so it can be matched with an ENST from the cds dictionary)
        mutation_ind_one = list((map(int,(re.findall("\d+",(row [17].partition("_")[0]))))))
        mutation_ind_two = list((map(int,(re.findall("\d+",(row [17].partition('_')[2]))))))
        mut_dict += [[row[4], gene_ide, mutation_ind_one,mutation_ind_two]]
for gn, gc in genedict.items():
    for s in mut_dict:
        if (gn == s[1]) and (len(gc) - s[2][0] > 0):                               #this loop runs if the mutation ENST and cds ENST match
            genecode = ''
            gene_mut_pos[gn] = s[2], s[3]

            for q in range(len(gc)):
                if len(s[3]) == 0:                  #creates a new string which takes values from the original cds as long as it is not at the mutation location
                    if q == s[2][0]:
                        genecode += ''                          #if it is at the mutation location then it returns a blank, i.e. that base is deleted
                    elif q != s[2][0]:
                        genecode += gc [q]
                elif len(s[3]) != 0:
                    if int(s[2][0]) < q < int(s[3][0]):           # this allows long sections of bases to be deleted from the cds
                        genecode += ''
                    elif int(s[2][0]) < q or int(s[3][0]) > q:
                        genecode += gc [q]
            rna = ''
#            for j in genecode:                                  # take the mutated cds and translate it to mRNA
#                if j == 'T':
#                    rna += 'U'
#                elif j != 'T':
#                    rna += j
            codon = ''
            for i in range(0, len(rna), 3):                     # loop through the rna string 3 items at a time
                codons = (rna[i: i + 3])                        # take the three items and store them as a variable 'codons'
                if codons in codon_list:                        # run codons through the codon dictionary and if they match to a key, add the corresponding value to 'codon'
                    codon += codon_list[codons]                 # this builds up a chain of single letter codes refelecting a peptide sequence in the variable 'codon'
                elif codons not in codon_list:
                    codon += ''
            if 'Z' in codon:                                    # find the position where the first stop codon occurs, the cds would be terminated at this point
                stopcod = codon.index('Z')
            if int(s[2][0])%3 == 0:
                prevnuc = [str(genecode [int(s[2][0]) -4 : int(s[2][0])-1])]
            elif s[2][0]%3 == 1:
                prevnuc = [str(genecode [int(s[2][0]) - 5 : int(s[2][0]) -2]) + str(genecode[int(s[2][0])-1])]
            elif s[2][0]%3 == 2:
                prevnuc = [str(genecode [int(s[2][0]) - 6 : int(s[2][0]) -3]) + str(genecode[int(s[2][0])-2]) + str(genecode[int(s[2][0])-1])]
            postnuc = str(genecode[int(s[2][0]) + 1 : int(s[2][0]) + 4])
            #perc = (stopcod / len(codon)) * 100
            diff = (stopcod * 3) - s[2][0]
            #aa_mut_pos[s[0]] = [gn, stopcod, prevnuc, len(codon), diff]      # make a dictionary so each ENST has a corresponding stopcodon
            #print(s[0], gn, stopcod, prevnuc, len(codon), diff)
            #output.write(postnuc + '\n')
            print(postnuc)
            s = []