# penncnv2bed.py
# C: Feb 10, 2016
# M: Dec 26, 2016
# A: Leandro Lima <lelimaufc@gmail.com>

# Program to create BED file from PennCNV output (.rawcnv) to be open on IGV or UCSC genome browser.
# Usage: python penncnv2bed.py result.rawcnv > result.bed


import sys

def main():

    print "track name=PennCNV-Seq"
    filename = sys.argv[1]

    lines = open(filename).read().split('\n')
    while lines[-1] == '':
        lines.pop()

    # chr1:12905142-13010536        numsnp=18     length=105,395     state1,cn=0 genome.baflrr startsnp=chr1:12905142-12905142 endsnp=chr1:13010536-13010536
    for line in lines:
        region, numsnp, length, state_cn, name, startsnp, endsnp = line.split()
        chrom  = region.split(':')[0]
        start  = region.split(':')[1].split('-')[0]
        end    = region.split(':')[1].split('-')[1]
        copies = state_cn.split(',cn=')[1]
        if int(copies) == 0:
            color = '128,0,0'
        elif int(copies) == 1:
            color = '255,0,0'
        elif int(copies) == 3:
            color = '0,0,255'
        elif int(copies) == 4:
            color = '0,0,128'
        elif int(copies) == 2:
            color = '128,128,128'
        print '%s\t%s\t%s\tcn=%s\t0\t*\t%s\t%s\t%s' % (chrom, start, end, copies, start, end, color)


if __name__ == '__main__':
    main()


