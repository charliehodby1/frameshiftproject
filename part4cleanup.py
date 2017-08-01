# make sense of delalldomdat.txt
import os
# os.chdir('/mnt/nfs2/biochem/hb280')
with open('delalldomdat.txt', 'r') as dadd:
    for line in dadd:
        line = line.split(',')
        ENST = line[0][4:-1]
        stopcod = line[2][3:-2]
        print(ENST, stopcod)
        domains = line[5:]
        b = 0
        for a in domains:
            b = b+1
            if b%3 != 0:
                if a.startswith(' [[['):
                    a = a[6:-1]
                    print(a)
                elif a.startswith(' u'):
                    a = a[3:-2]
                    print(a)
                elif a.startswith(' [[u'):
                    a = a[5: -1]
                    print(a)
            elif b%3 == 0:
                print(a[3:-4])