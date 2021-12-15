#!/usr/bin/env python

def filter_sam(inFile):
	"""
	Filter unaligned reads, reads with secondary alignments, reads with MAPQ < 30
	"""
	samfile = pysam.AlignmentFile(inFile, "r")

	primary = {}
	secondary = {}

	for read in samfile.fetch(until_eof=True):
		if not read.is_unmapped: # Read mapped
			if read.is_secondary or read.is_supplementary: # Read secondary
				if read.qname in primary.keys(): # Read has primary alignment
					if read.reference_length >= (0.8 * primary[read.qname].reference_length): # Secondary alignment length is longer than 80% primary alignemnt length
						secondary[read.qname] = primary[read.qname] # Move primary to secondary
						del primary[read.qname] # Remove primary
			else: # Read primary
				if read.mapping_quality > 30: # Read good mapping quality
					primary[read.qname] = read 
				else: # Read bad mapping quality
					secondary[read.qname] = read
		else:
			secondary[read.qname] = read # Read unmapped

	return [primary, secondary] 

def output_primary(primary, outF, inF):
	"""
	Output primary alignments as sam.
	"""
	samfile = pysam.AlignmentFile(inF, "r")
	pFile = pysam.AlignmentFile(outF, "w", template=samfile)

	for v in primary.values():
		pFile.write(v)

def output_secondary(secondary, outF):
	"""
	Output reads with secondary alignments, unaliged or aligned with MAPQ < 30 as fasta.
	"""
	sFile = open(outF, "w")

	for v in secondary.values():
		sFile.write(">{}\n{}\n".format(v.qname[:50], v.query_sequence))


if __name__ == "__main__":
	import pysam
	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument("-i", help="Input sam file.", action="store", required=True)
	parser.add_argument("-p", help="Output sam file name for primary alignments.", action="store", required=True)
	parser.add_argument("-s", help="Output fasta file name for secondary alignments.", action="store", required=True)
	args = parser.parse_args()

	primary, secondary = filter_sam(args.i)
	output_primary(primary, args.p, args.i)
	output_secondary(secondary, args.s)
