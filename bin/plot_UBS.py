#!/usr/bin/env python

def parse_fasta(fasta):
	"""
	Parse FASTA file and retrieve ids
	"""
	ids = []

	for seq in SeqIO.parse(fasta, "fasta"):
		ids.append(len(seq.id))

	return ids 

def classify_sizes(ids, names, sizes):
	"""
	Classify sizes by ids
	"""
	selected = []
	not_selected = []

	for n, s in zip(names, sizes):
		if n in ids:
			selected.append(s)
		else:
			not_selected.append(s)

	return (selected, not_selected)

def plot(sel, no_sel, sample):
	"""
	Plot an interactive histogram of umi bin sizes
	https://plotly.com/python/histograms/
	"""
	x0 = sel
	x1 = no_sel

	fig = go.Figure()
	fig.add_trace(go.Histogram(x=x0, name='filtered & trimmed UMIs', marker_color='#70a5fa'))
	fig.add_trace(go.Histogram(x=x1, name='discarded UMIs', marker_color='#3f4c61'))

	# The two histograms are drawn on top of another
	fig.update_layout(barmode='stack')
	fig.write_html(sample + "_ubs.html")


if __name__ == "__main__":
	from Bio import SeqIO
	import plotly.graph_objects as go
	import argparse
	import random
	import numpy

	parser = argparse.ArgumentParser()
	parser.add_argument("-s", help="Sample name", action="store", required=True)
	parser.add_argument("-t", help="Fasta file with trimmed reads", action="store", required=True)
	parser.add_argument("--sizes", help="UMI cluster sizes", nargs="*", action="store", type=int, required=True)
	parser.add_argument("--names", help="UMI cluster names", nargs="*", action="store", type=str, required=True)
	args = parser.parse_args()

	trimmed = parse_fasta(args.t)
	sel, no_sel = classify_sizes(trimmed, args.names, args.sizes)
	plot(sel, no_sel, args.s)
