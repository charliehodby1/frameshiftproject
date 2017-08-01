# the one i actually used adding dom data to domrecdeloutput
import os
import csv
import re
import gzip
from collections import defaultdict
# os.chdir("/mnt/nfs2/biochem/hb280")
enst_mut_stop = []
with open('domrecdeloutput.txt', 'r') as drdo:
    for line in drdo:
        line = line.split(' " ')
        line = line[0]
        ENST = (line[3 : 18])
        line = line.split('"')
        if len(line) > 1:
            mutname = (line[1][1:])
            stopcod = (line[2][3: -2])
        if len(ENST) > 5:
            enst_mut_stop += [[ENST, mutname, stopcod]]

domain_dict = {}                                                    # open a dictionary ('domain_dict') to store results from the ENST-Pfam reference file
enpf = open('enstgenenamenew.txt','r')
for entry in csv.reader(enpf, dialect = 'excel-tab'):                                        # this means we dont consider the first piece of data as it is header information
        tremblacc = entry[3]                                   # pull out the Pfam reference
        swiss = entry[4]
        ENST = entry[1]                                         # pull out the ENST reference
        domain_dict[ENST] = swiss,tremblacc                    # combine the data so that each ENST has a corresponding Pfam reference in dictionary 'domain_dict'
dom_dat = {}                                                      # open a second dictionary which will eventually contain an ENST reference corresponding to domain data
pfam = gzip.open('9606.tsv.gz', 'rb')
for info in pfam:
    info = info.decode('utf-8')                                     # decode the data from compressed file
    info = info.strip().split(' ')                                  # split the data so we can more easily consider individual elements
    pflist = []
    for d in info:
        pflist.append(d.split('\t'))                                # remove the '\t' from the data and store each line in a list
    newlist = pflist[0]                                             # extract the lines of data from the list
    if len(newlist) > 1 and newlist[7] == 'Domain':                 # exclude the first line as this contains header infromation
        pfamref = newlist[0]                                        # identify the Pfam reference
        domain_bound = [newlist[1], newlist[2]]                     # identify the domain boundaries
        domain_name = newlist[6]                                    # identify the domain name
        dom_dat[pfamref] = dom_dat.get(pfamref, []) + [[domain_bound, domain_name]]# create a dictionary where an ENST reference corresponds to domain information
ems_domdict = []


for enst, dom in domain_dict.items():
    for a in enst_mut_stop:
        if a[0] == enst:
            ems_domdict += [[a, dom]]
            a = []

k = 0            
ems_dict_dom = []
for pfamr, domdata in dom_dat.items():
    for b in ems_domdict:
        if b[1][0] == pfamr:
            k = k + 1
            ems_dict_dom += [[b, domdata]]
            print(b, domdata)
            b = []
