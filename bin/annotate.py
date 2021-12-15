#!/usr/bin/env python

def collect_peaks(bedFile):
	"""Collect all peaks from bed file and look for repeats"""
	# Inititalyze Dfam API
	url = "https://dfam.org/api/annotations"
	params = {
		"assembly": "hg38",
		"nrph": "true",
	}
	# Collect peaks
	peak_file = open(bedFile)
	annotations = {}
	for line in peak_file:
		line = line.strip()
		line = line.split("\t")
		params["chrom"] = line[0]
		params["start"] = int(line[1]) - 50
		params["end"] = int(line[2]) + 50
		# Retrieve annotations
		response = requests.get(url, params=params)
		annotations[response.json()["query"]] = {}
		annotations[response.json()["query"]]["hits"] = response.json()["hits"]
		annotations[response.json()["query"]]["tandem_repeats"] = response.json()["tandem_repeats"]

	return(annotations)

def get_counts(string, n):
	string = string.split()
	return string[n]

def parse_repeatmasker(tbl):
	labels=["ALUs", "MIRs", "LINE1", "LINE2", "L3/CR1", "ERVL", "ERVL-MaLRs", "ERV_classI", "ERV_classII", "hAT-Charlie", "TcMar-Tigger", "", "", "", "", ""]
	parents=["SINEs", "SINEs", "LINEs", "LINEs", "LINEs", "LTR elements", "LTR elements", "LTR elements", "LTR elements", "DNA elements", "DNA elements", "Unclassified", "Small RNA", "Satellites", "Simple Repeats", "Low complexity"]
	counts = []
	file = open(tbl)

	for line in file:
		if line.startswith("SINEs"):
			for i in range(2):
				line = file.readline()
				counts.append(get_counts(line, 1))
			for i in range(2):
				line = file.readline()
			for i in range(3):
				line = file.readline()
				counts.append(get_counts(line, 1))
			for i in range(2):
				line = file.readline()
			for i in range(4):
				line = file.readline()
				counts.append(get_counts(line, 1))
			for i in range(2):
				line = file.readline()
			for i in range(2):
				line = file.readline()
				counts.append(get_counts(line, 1))
			line = file.readline()
			line = file.readline()
			counts.append(get_counts(line, 1))
			for i in range(4):
				line = file.readline()
			line = file.readline()
			counts.append(get_counts(line, 2))
			line = file.readline()
			line = file.readline()
			counts.append(get_counts(line, 1))
			for i in range(2):
				line = file.readline()
				counts.append(get_counts(line, 2))

	df = pd.DataFrame(dict(type=parents, query=labels, counts=counts))

	return df


def parse_annotations(annotations, post):
	"""Retrieve and format annotations"""
	queries = []
	types = []
	tandem = []
	notrep = 0
	rep = 0
	for k in annotations.keys():
		if len(annotations[k]["hits"]) >= 1:
			queries.append([h["query"] for h in annotations[k]["hits"]]) 
			types.append([h["type"] for h in annotations[k]["hits"]]) 
		if len(annotations[k]["tandem_repeats"]) >= 1:
			tandem.append([h["type"] for h in annotations[k]["tandem_repeats"]]) 
		if len(annotations[k]["hits"]) < 1 and len(annotations[k]["tandem_repeats"]) < 1:
			if post:
				queries.append(["Non repetitive"]) 
				types.append([""]) 
			notrep += 1
		else:
			rep += 1
		
	# append tandem repeats to hits
	[queries.append(t) for t in tandem]
	types.append(["tandem_repeat"]*len(tandem))
	# Format list to remove sublists
	queries = [item for sublist in queries for item in sublist]
	types = [item for sublist in types for item in sublist]
	# Count annotations
	joint = [(q, t) for q, t in zip(queries, types)]
	counts = [joint.count(e) for e in set(joint)]
	# Format annotations to data frame to plot
	df = pd.DataFrame(dict(query=[s[0] for s in set(joint)], type=[s[1] for s in set(joint)], counts=counts))
	r_count = [["Repetitive", "Non repetitive"], [rep, notrep]]

	return([df, r_count])

def plot(df, r_count, outFile, post):
	"""Plot regions to html pie chart"""
	if not df.empty:
		sunburst = px.sunburst(df, path=['type', 'query'], values='counts', 
			color='type', color_discrete_sequence=px.colors.qualitative.Set3)
		if not post:
			fig = make_subplots(rows=1, cols=2, specs=[[{'type':'domain'}, {'type':'domain'}]])
			fig.add_trace(go.Pie(values=r_count[1], labels=r_count[0]), 1, 1)
			fig.add_trace(go.Sunburst(labels=sunburst['data'][0]['labels'].tolist(),
				parents=sunburst['data'][0]['parents'].tolist(),
				marker=dict(
					colors=px.colors.qualitative.Set3)),
				1, 2)
			fig.update_layout(uniformtext=dict(minsize=20, mode='show'))
			fig.write_html(outFile)
		else:
			sunburst.update_layout(uniformtext=dict(minsize=20, mode='show'))
			sunburst.write_html(outFile)


if __name__ == "__main__":
	from plotly.subplots import make_subplots
	import plotly.graph_objects as go
	import plotly.express as px
	import pandas as pd
	import requests
	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument("-i", help="Input bed/tbl file", action="store", required=True)
	parser.add_argument("-o", help="Output html file name", action="store", required=True)
	parser.add_argument("--post", help="Bed file is post peak calling", action="store_true")
	parser.add_argument("--tbl", help="Input file is a tbl file, RepeatMasker output", action="store_true")
	args = parser.parse_args()

	if args.tbl:
		df = parse_repeatmasker(args.i)
		plot(df, "", args.o, True)
	else:
		dfs = parse_annotations(collect_peaks(args.i), args.post)
		plot(dfs[0], dfs[1], args.o, args.post)

