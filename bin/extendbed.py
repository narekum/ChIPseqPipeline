# Author: pzs
# Extended by McBryan to add chrmEnd trimming

import csv
from genemapping.chrmEnds import ChromosomeEnds

def extendOne(fname, extendlen):
	reader = csv.reader(open(fname), delimiter="\t")
	outfile = fname + ".extended"
	print "writing to outfile", outfile
	writer = csv.writer(open(outfile, "w"), delimiter="\t")
	for row in reader:
		chr = row[0]
		start = int(row[1])
		end = int(row[2])
		strand = row[5]
		if strand == "+":
			row[2] = min(ends[chr],start + extendlen)
		elif strand == "-":
			row[1] = max(0,end - extendlen)
		
		writer.writerow(row)

if __name__ == "__main__":
	from glob import glob
	import sys
	if len(sys.argv) != 4:
		print "usage: extendbed.py <beddir> <genomebuild> <extendlength>"
		sys.exit(1)
		
	ends = ChromosomeEnds(sys.argv[2])
	
	beddir = sys.argv[1]
	extendlen = int(sys.argv[3])
	allfiles = glob(beddir + "/*.bed")
	for fname in allfiles:
		extendOne(fname, extendlen)
