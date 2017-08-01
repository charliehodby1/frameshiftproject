
import os
# os.chdir("/mnt/nfs2/biochem/hb280")
delo_li = []
dr_li = []
with open ('indoutput.txt', 'r') as delo:
    for entry in delo:
        entry = entry.split(',')
        markerdelo = entry[1][2 :17]
        markervalue = entry[0], entry[2], entry[3], entry[4]
        delo_li += [[markerdelo, markervalue]]


with open('domrec.txt', 'r') as dr:
    for line in dr:
        line = line.split(',')
        markerdr = line[0]
        valuedr = line[1]
        dr_li += [[markerdr, valuedr]]

for b in dr_li:
    for a in delo_li:
        if a[0] == b[0]:
            print((str(b[1][: 2]) + ',' + str(a[0]) + ',' +  str(a[1])) + '\n')
