"""
gsToBed.py

Convert GUIDE-Seq to BED file for easy use with bedtools

"""
import sys

def gsToBed(inpath, outpath):
    with open(inpath, 'r') as infile:
        with open(outpath, 'w') as outfile:
            infile.readline()
            for line in infile:
                line = line.split('\t')
                if float(line[12]) > 0:
                    coordinates = line[0].replace(':', '-').split('-')
                    chromosome = coordinates[0]
                    start = coordinates[1]
                    end = coordinates[2]
                    name = line[3]
                    outfile.write('\t'.join([chromosome, start, end, name]))
                    outfile.write('\n')

if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    gsToBed(infile, outfile)