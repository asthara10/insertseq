#!/usr/bin/env python3

def parseInput(inFile):
	"""
	Parse RepeatMasker output to retrieve teh best hit (Repeat) for each read (Query) 
	which includes the start or end of the read.
	"""
	file = open(inFile)
	next(file) # Skip first 3 lines
	next(file)
	next(file)
	qname = "" # Query name == Read name
	rname = "" # Repeat name
	score = "" # Match score
	hits = {}
	
	for line in file:
		line = line.strip().split()
		#if int(line[5]) <= 1 or int(line[7][1:-1]) == 0: # Query match start <= 1 or num bases after mathc == 0
		if int(line[5]) <= 20 or int(line[7][1:-1]) <= 20: # Query match start <= 20 (repeat starts before the first 20 bp of the read) or num bases after mathc == 0 (repeat ends after last 20 bp of the read)
			if line[4] in hits.keys(): # Existent query
				if line[9] != hits[line[4]][0]: # Repeat is not the same as the one already assigned to query
					if int(line[0]) > hits[line[4]][1]: # Repeat score is higher than the already assigned
						hits[line[4]] = [line[9], int(line[0])] # Substitute assigend repeat
			else: # New read
				hits[line[4]] = [line[9], int(line[0])] # Add query name as key, a list of repeat name and alignment score as value	
	
	return hits

def countRepeats(hits):
	"""
	Count the occurences of each Repeat.
	Retrieve reads containing the repeat.
	"""
	counts = {}
	reads = {}
	
	for k, v in hits.items():
		try:
			counts[v[0]] += 1
			reads[v[0]].append(k)
		except KeyError:
			counts[v[0]] = 1
			reads[v[0]] = [k]
	
	counts_sorted = {k: v for k, v in sorted(counts.items(), key=lambda item: -item[1])}
	
	return [counts_sorted, reads]

def alignRepeats(reads, fasta):
	"""
	Perform a multiple sequence alignment for read in each repeat.
	"""
	reads_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
	for k, v in reads.items():
		if len(v) > 2:
			reads_list = []
			for r in v:
				reads_list.append(reads_dict[r])
			k = ''.join(e for e in k if e.isalnum())
			SeqIO.write(reads_list, k+"_temp.fasta", "fasta")
			cmd = ClustalwCommandline("clustalw2", infile=k+"_temp.fasta")
			stdout, stderr = cmd()
			align = AlignIO.read(k+"_temp.aln", "clustal")
			print(align)


def significantRepeats(repeats):
	"""
	Return only significant repeats.
	The ones that may contain an insertion.
	"""


def outputAllRepeatCounts(repeats, outf):
	"""
	Print to output file all found Repeats with count of reads containing that repeat
	"""
	out = open(outf, 'w')
	out.write("Repeat Name\tRead Count\n")
	
	for k, v in repeats.items():
		out.write("{}\t{}\n".format(k,v))



if __name__ == "__main__":
	from Bio.Align.Applications import ClustalwCommandline
	from plotly.subplots import make_subplots
	from Bio import SeqIO
	from Bio import AlignIO
	import plotly.graph_objects as go
	import plotly.express as px
	import pandas as pd
	import requests
	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument("-i", help="Input annotation file from RepeatMasker output", action="store", required=True)
	parser.add_argument("-f", help="Input fasta file with unmapped repeats", action="store", required=True)
	parser.add_argument("-o", help="Output called unanchored peaks file name", action="store", required=True)
	args = parser.parse_args()

	repeats, reads = countRepeats(parseInput(args.i))
	alignRepeats(reads, args.f)
	#outputAllRepeatCounts(repeats, args.o)
