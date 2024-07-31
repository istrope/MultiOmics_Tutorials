import re
import sys
import random
import string

def main():

	INFILE = open("hg38.goodzones.bwa.k50.bed", "r")
	OUTFILE = open("hg38.chrom.mappable.bwa.k50.txt", "w")

	prevChrom = "chr1"
	prevChromMappable = 0

	for x in INFILE:
		arow = x.rstrip().split("\t")
		thisChrom = arow[0]
		thisLength = int(arow[2]) - int(arow[1])

		if thisChrom == prevChrom:
			prevChromMappable += thisLength
		else:
			OUTFILE.write(prevChrom)
			OUTFILE.write("\t")
			OUTFILE.write(str(prevChromMappable))
			OUTFILE.write("\n")
			prevChromMappable = thisLength
		prevChrom = thisChrom

	OUTFILE.write(prevChrom)
	OUTFILE.write("\t")
	OUTFILE.write(str(prevChromMappable))
	OUTFILE.write("\n")
	
	INFILE.close()	
	OUTFILE.close()


if __name__ == "__main__":
	main()
