#!/usr/bin/env python3

def filter_repeats(infile, threshold):
	"""Filter all repeats that pass the coverage threshold"""
	tbl = open(infile)

	repeats = {}
	filtered = {}

	for line in tbl:
		if not line.startswith("#"):
			line = line.strip().split()
			try:
				repeats[line[0]]["reads"].append(line[2])
			except KeyError:
				repeats[line[0]] = {}
				repeats[line[0]]["reads"] = [line[2]]
				repeats[line[0]]["accession"] = [line[1]]
				repeats[line[0]]["description"] = [line[-1]]

	for rep in repeats.keys():
		if len(repeats[rep]["reads"]) >= threshold and len(repeats[rep]["reads"]) > 2:
			filtered[rep] = repeats[rep]

	return filtered


def obtain_reads(repeats, fasta):
	"""obtain repeats read sequences from original fasta file"""
	reads = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

	for rep in repeats.keys():
		sequences = ""

		for read in repeats[rep]["reads"]:
			sequences += reads[read].format("fasta")

		with open("tmp/" + rep + "_sequences.fasta", "w") as output:
			output.write(sequences)


def cluster(target):
	file = "tmp/" + target + "_sequences.fasta"
	command = ["vsearch", "--cluster_size", file, "--clusters", "tmp/" + target + "_cluster", "--id", "0.8", "--clusterout_id"]
	p = subprocess.Popen(command, stdout=subprocess.PIPE)
	output, err = p.communicate()


def filter_clusters(target, threshold):
	clusters = []

	cluster_files = [join("tmp", f) for f in listdir("tmp") if f.startswith(target) and "." not in f and isfile(join("tmp", f))]

	for cluster in cluster_files:
		cluster_reads = list(SeqIO.parse(cluster, "fasta"))

		if len(cluster_reads) >= threshold:
			new_reads = []

			for r in cluster_reads:
				new_reads.append(r.id)

			print("Cluster {} has {} reads and PASS the coverage threshold of {}.".format(cluster, len(cluster_reads), threshold))
			clusters.append(new_reads)
		
		else:
			print("Cluster {} has {} reads and does NOT PASS the coverage threshold of {}.".format(cluster, len(cluster_reads), threshold))

	for filename in os.listdir("tmp"):
		#if not filename.endswith(".fasta"):
		if "." not in filename:
			file_path = os.path.join("tmp", filename)
			os.unlink(file_path)

	return clusters


def perform_msa(target, cluster, read_names, fasta):
	"""Obtain multiple sequence alignment for each target cluster"""
	reads = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
	sequences = ""

	for read in read_names:
		sequences += reads[read].format("fasta")

	fast_output = "tmp/" + target + "_cluster" + str(cluster) + "_reads.fasta"
	with open(fast_output, "w") as output:
		output.write(sequences)

	clust_output = "tmp/" + target + "_cluster" + str(cluster) + "_msa.fasta"
	command = ["clustalo", "-i", fast_output, "-o", clust_output]
	p = subprocess.Popen(command, stdout=subprocess.PIPE)
	output, err = p.communicate()


"""def perform_alignment(target, cluster, read_names, fasta, ref_path):
	#Align each cluster to it's target consensus reference
	reads = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
	reference = ref_path + target + ".fasta"
	sample = "tmp/" + target + "_cluster" + str(cluster)
	sequences = ""
	sam = sample+".sam"
	bam = sample+".sorted.bam"
	bai = sample+".sorted.bam.bai"
	msa = sample+".msa"

	for read in read_names:
		sequences += reads[read].format("fasta")

	fast_output = sample + "_reads.fasta"
	with open(fast_output, "w") as output:
		output.write(sequences)

	command = ["minimap2", "-t", "4", "-ax", "map-ont", reference, fast_output, "-o", sam]
	p = subprocess.Popen(command, stdout=subprocess.PIPE)
	output, err = p.communicate()
	command = ["samtools", "sort", "-@", "4", "-O", "BAM", "-o", bam, sam]
	p = subprocess.Popen(command, stdout=subprocess.PIPE)
	output, err = p.communicate()
	command = ["samtools", "index", "-@", "4", bam, bai]
	p = subprocess.Popen(command, stdout=subprocess.PIPE)
	output, err = p.communicate()
	command = ["bam2msa", "--display-left-softclip", "--display-right-softclip", reference, bam]
	p = subprocess.Popen(command, stdout=subprocess.PIPE)
	output, err = p.communicate()
	with open(msa, 'w') as out_msa:
		out_msa.write(str(output))
"""

def filter_sds(target, cluster, tgen, tdif):
	"""Obtain standard deviations from msa and filter"""
	msa = "tmp/" + target + "_cluster" + str(cluster) + "_msa.fasta"
	#msa = "tmp/" + target + "_cluster" + str(cluster) + ".msa"
	alignment = SeqIO.parse(msa, "fasta")
	pos_start = []
	pos_end = []

	for read in alignment:
		pos_start.append(len(read.seq) - len(read.seq.lstrip('-')))
		pos_end.append(len(read.seq) - len(read.seq.rstrip('-')))

	sd_start = np.std(pos_start)
	sd_end = np.std(pos_end)
	if sd_start > sd_end:
		sd_gen = sd_start
		sd_ins = sd_end
	else:
		sd_gen = sd_end
		sd_ins = sd_start

	if sd_gen > tgen and (sd_gen - sd_ins) > tdif:
		return True
	else:
		return False


def filter_target(target, threshold, fasta, filt_clusters, tgen=20.54449, tdif=17.66657):
	cluster(target)
	clusters = filter_clusters(target, threshold)

	if len(clusters) > 0 and len(clusters[0]) > 0:
		print("Target {} has {} clusters, and does PASS to next filtering steps.".format(target, len(clusters)))
		if len(clusters) == 1:
			perform_msa(target, 0, clusters[0], fasta)
			#perform_alignment(target, 0, clusters[0], fasta, "/data/ReferenceGenomes/hg38/repeats/")
			if filter_sds(target, 0, tgen, tdif):
				filt_clusters[target + "_cluster0"] = clusters[0]
		else:
			for num, clust in enumerate(clusters):
				perform_msa(target, num, clust, fasta)
				#perform_alignment(target, num, clust, fasta, "/data/ReferenceGenomes/hg38/repeats/")
				if filter_sds(target, num, tgen, tdif):
					filt_clusters[target + "_cluster" + str(num)] = clust
		return filt_clusters
	else:
		print("Target {} has {} clusters, and does NOT PASS to next filtering steps.".format(target, len(clusters)))
		return filt_clusters


def print_output(peaks, output):
	"""Print bed file containing only insertions passing shape filter"""
	with open(output, 'w') as out:
		for repeat in peaks.keys():
			out.write("{}\t{}\t{}\n".format(repeat, len(peaks[repeat]), '\t'.join(peaks[repeat])))


if __name__ == "__main__":
	import pysam
	import numpy as np
	import argparse
	import os
	import subprocess
	import numpy as np
	from Bio import pairwise2
	from Bio import SeqIO
	from Bio.Seq import Seq
	from os import listdir
	from os.path import isfile, join

	parser = argparse.ArgumentParser()
	parser.add_argument("--fasta", help="Fasta file with unanchored reads", action="store", required=True)
	parser.add_argument("--tbl", help="hmmscan output tbl file", action="store", required=True)
	parser.add_argument("--tcov", help="Minimum number of reads per peak", action="store", required=True, type=int)
	parser.add_argument("-o", help="Output text file with repeats passing peak filter", action="store", required=True)
	parser.add_argument("--tdif", help="Threshold of standard deviation differences for peak selection. [default = 207.5480]", action="store", type=float)
	parser.add_argument("--tgen", help="Genome standard deviation threshold for peak selection. [default = 222.7823]", action="store", type=float)    
	args = parser.parse_args()

	repeats = filter_repeats(args.tbl, int(args.tcov))
	obtain_reads(repeats, args.fasta)
	filt_clusters = {}

	for target in repeats.keys():
		if args.tgen is None and args.tdif is None:
			filt_clusters = filter_target(target, int(args.tcov), args.fasta, filt_clusters)
		elif args.tgen is None:
			filt_clusters = filter_target(target, int(args.tcov), args.fasta, filt_clusters, tdif=args.tdif)
		elif args.tdif is None:
			filt_clusters = filter_target(target, int(args.tcov), args.fasta, filt_clusters, tgen=args.tgen)
		else:
			filt_clusters = filter_target(target, int(args.tcov), args.fasta, filt_clusters, args.tgen, args.tdif)

	print_output(filt_clusters, args.o)
