#!/usr/bin/env python

def get_reference_vector():
    """Obtain reference vector from binary file.
    Reference vector has beed calculated withthe mean per position of 42 true insertions
    selected manually."""
    reference_vector = np.load("reference_vector.npy")
    return(reference_vector)

def parse_bed(bedFile, LOD=2):
    """Obtain peak regions from bed file.
    Filter insertions with enough coverage.
    LOD: The required coverage for an insertion to adjust the expected distribution."""
    insertions = {}
    with open(bedFile) as bed:
        for line in bed:
            line = line.strip()
            line = line.split("\t")
            if int(line[3]) > LOD: # Insertions with more than LOD coverage
                insertions[line[0]+":"+line[1]+"-"+line[2]] = [line[0], int(line[1]), int(line[2]), int(line[3])] # chromosome, start, end, coverage
    return(insertions)

def obtain_start_end(bam, insertions):
    """Obtain start and end positions of all reads from each insertion."""
    positions = {}
    bamFile = pysam.AlignmentFile(bam, "rb")
    for ins in insertions.keys():
        reads = bamFile.fetch(insertions[ins][0], insertions[ins][1], insertions[ins][2], until_eof=True)
        positions[ins] = {"start":[], "end":[]}
        for r in reads:
            positions[ins]["start"].append(r.reference_start)
            positions[ins]["end"].append(r.reference_end)
        positions[ins]["start"] = np.array(positions[ins]["start"])
        positions[ins]["end"] = np.array(positions[ins]["end"])
    return(positions)

def reject_outliers(data, m=2):
    """Reject """
    return data[abs(data - np.mean(data)) < m * np.std(data)]

def check_positions(positions, insertions, threshold=100):
    """Check that all reads form an insertion start at the same positions and end in differnt positions, or viceversa."""
    selected = []
    for ins in insertions.keys():
        if ((np.std(reject_outliers(positions[ins]["start"])) < threshold and 
            np.std(reject_outliers(positions[ins]["end"])) > threshold) or
                (np.std(reject_outliers(positions[ins]["start"])) > threshold and
                np.std(reject_outliers(positions[ins]["end"])) < threshold)):
            selected.append(ins)
    return(selected)

def check_positions_outlierfilter(positions, insertions, threshold=100):
    """Check that all reads form an insertion start at the same positions and end in differnt positions, or viceversa."""
    selected = []
    for ins in insertions.keys():
        if ((np.std(reject_outliers(positions[ins]["start"])) < threshold and 
            np.std(reject_outliers(positions[ins]["end"])) > threshold) or
                (np.std(reject_outliers(positions[ins]["start"])) > threshold and
                np.std(reject_outliers(positions[ins]["end"])) < threshold)):
            selected.append(ins)
    return(selected)

def check_positions_outlierfilter_twovals(positions, insertions, threshold_genome=178, threshold_insertion=23):
    """Check that all reads form an insertion start at the same positions and end in differnt positions, or viceversa."""
    selected = []
    for ins in insertions.keys():
        if ((np.std(reject_outliers(positions[ins]["start"])) < threshold_genome and 
            np.std(reject_outliers(positions[ins]["end"])) > threshold_insertion) or
                (np.std(reject_outliers(positions[ins]["start"])) > threshold_genome and
                np.std(reject_outliers(positions[ins]["end"])) < threshold_insertion)):
            selected.append(ins)
    return(selected)

def check_positions_twovals(positions, insertions, threshold_genome=236, threshold_insertion=202):
    """Check that all reads form an insertion start at the same positions and end in differnt positions, or viceversa.
    threshold_genome is the minimum standard deviation of selected true insertions SD at genome site.
    threshold_insertion is the maximum SD of selected true insertions SD at insertion site."""
    selected = []
    orientation = []
    for ins in insertions.keys():
        if (np.std(positions[ins]["start"]) < threshold_genome and np.std(positions[ins]["end"]) > threshold_insertion):
            selected.append(ins)
            orientation.append(1) # Insertion is happening at start
        elif (np.std(positions[ins]["start"]) > threshold_genome and np.std(positions[ins]["end"]) < threshold_insertion):
            selected.append(ins)
            orientation.append(2) # Insertion is happening at end
    return((selected, orientation))

def check_positions_difference(positions, insertions, threshold=207.54804):
    """Check that all reads form an insertion start at the same positions and end in differnt positions, or viceversa.
    Using the difference between genome-insertion standard deviations."""
    selected = []
    orientation = []
    for ins in insertions.keys():
        if abs( np.std(positions[ins]["start"]) - np.std(positions[ins]["end"]) ):
            selected.append(ins)
            if np.std(positions[ins]["start"]) < np.std(positions[ins]["end"]):
                orientation.append(1) # Insertion is happening at start
            else:
                orientation.append(2)
    return((selected, orientation))

def check_positions_combination(positions, insertions, threshold_genome=222.7823, threshold_difference=207.5480):
    """Check that all read from an insertion start at the same positions and end in different positions, or viceversa.
    Use threshold of standard deviation at genome side and difference between SD at genome and insertion side."""
    selected = []
    orientation = []
    for ins in insertions.keys():
        sd_st = np.std(positions[ins]["start"]) # Standard deviation of read starts
        sd_end = np.std(positions[ins]["end"]) # Standard deviation of read ends
        if sd_st < sd_end: # start == insertion side & end == genome side
            sd_ins = sd_st
            sd_gen = sd_end
            orient = 1  # Insertion is happening at start
        else:
            sd_ins = sd_end
            sd_gen = sd_st
            orient = 2  # Insertion is happening at end
        if (sd_gen > threshold_genome) & ((sd_gen - sd_ins) > threshold_difference):
            selected.append(ins)
            orientation.append(orient)
    return((selected, orientation))

def obtain_coverage_distribution(bam, insertions, selected):
    """Obtain coverage distribution per insertion"""
    coverage = {}
    bamFile = pysam.AlignmentFile(bam, "rb")
    for ins in selected:
        pileup = bamFile.pileup(insertions[ins][0], insertions[ins][1], insertions[ins][2], max_depth=1000000)
        coverage[ins] = [col.nsegments for col in pileup]
        coverage[ins] = np.array(coverage[ins])
    return(coverage)

def scale_distributions(coverage):
    """Scale coverage distributions from 0 to 1."""
    scaled = {}
    for ins in coverage.keys():
        scaled[ins] = coverage[ins]/max(coverage[ins])
    return(scaled)

def calculate_euclidean(scaled):
    """Calculate euclidean distace between a reference true vector of coverage 
    and the scaled coverage of each found insertion.
    Reference vector has beed calculated withthe mean per position of 42 true insertions
    selected manually."""
    reference = get_reference_vector()
    distances = {}
    for ins, array in scaled.items():
        if array.shape[0] > reference.shape[0]:
            reshaped = np.zeros(array.shape)
            reshaped[:reference.shape[0],] = reference
            distances[ins] = np.linalg.norm(array - reshaped)
        else:
            reshaped = np.zeros(reference.shape)
            reshaped[:array.shape[0],] = array
            distances[ins] = np.linalg.norm(reshaped - reference)
    return(distances)

def calculate_DHD(scaled):
    """Calculate Directed Hausdorff distace between a reference true vector of coverage 
    and the scaled coverage of each found insertion.
    Reference vector has beed calculated withthe mean per position of 42 true insertions
    selected manually."""
    reference = get_reference_vector()
    distances = {}
    for ins, array in scaled.items():
        if array.shape[0] > reference.shape[0]:
            reshaped = np.zeros(array.shape)
            reshaped[:reference.shape[0],] = reference
            u = np.reshape(array, (1, -1))
            v = np.reshape(reshaped, (1, -1))
        else:
            reshaped = np.zeros(reference.shape)
            reshaped[:array.shape[0],] = array
            u = np.reshape(reshaped, (1, -1))
            v = np.reshape(reference, (1, -1))
        distances[ins] = scipy.spatial.distance.directed_hausdorff(u, v)[0]
    return(distances)

def print_best_dist(coverage):
    best = []
    for ins in coverage.keys():
        dist = distfit()
        dist.fit_transform(coverage[ins])
        best.append(dist.summary.distr[0])
    return(best)

def shape_filter(coverage, distribution="gamma", threshold=0.0007547906866842447):
    """Select only insertions with propper coverage distribution"""
    filtered = []
    for ins in coverage.keys():
        fit = distfit(distr = distribution)
        fit.fit_transform(coverage[ins])
        if float(fit.summary.RSS) <= threshold:
            filtered.append(ins)
    return(filtered)

def calculate_inssite(coverage, insertions, filtered, orientation, outSite):
    """Calculate insertion happening site"""
    out = open(outSite, 'w')
    for ins, orient in zip(filtered, orientation):
        out.write("{}\t{}\t{}\n".format(insertions[ins][0], insertions[ins][orient], insertions[ins][orient]))

def print_output(insertions, filtered, output):
    """Print bed file containing only insertions passing shape filter"""
    with open(output, 'w') as out:
        for f in filtered:
            out.write("{0}\t{1}\t{2}\t{3}\n".format(*insertions[f]))

def print_output_distance(insertions, selected, distances, output):
    """Print bed file containing only insertions passing shape filter"""
    with open(output, 'w') as out:
        for true in selected:
            out.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(*insertions[true], distances[true])) # chromosome, start, end, coverage, distance


if __name__ == "__main__":
    import pysam
    from distfit import distfit
    import numpy as np
    import argparse
    import scipy
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", help="Bam file of mapped reads", action="store", required=True)
    parser.add_argument("--bed", help="Bed file with stignificant called peaks", action="store", required=True)
    parser.add_argument("--oFilt", help="Output bed file with peaks passing shape filter", action="store", required=True)
    parser.add_argument("--oSite", help="Output bed file with peaks passing shape filter and only insertion site/position", action="store", required=True)
    parser.add_argument("--tdif", help="Threshold of standard deviation differences for peak selection. [default = 207.5480]", action="store", type=float)
    parser.add_argument("--tgen", help="Genome standard deviation threshold for peak selection. [default = 222.7823]", action="store", type=float)
    #parser.add_argument("--tins", help="Insertion standard deviation threshold for peak selection. [default = 202]", action="store", type=float)
    args = parser.parse_args()

    insertions = parse_bed(args.bed, 10)
    positions = obtain_start_end(args.bam, insertions)

    #if args.t1 is None and args.t2 is None:
        #selected, orientation = check_positions_twovals(positions, insertions)
    #elif args.t1 is None:
        #selected, orientation = check_positions_twovals(positions, insertions, args.t2)
    #elif args.t2 is None:
        #selected, orientation = check_positions_twovals(positions, insertions, args.t1)
    #else:
        #selected, orientation = check_positions_twovals(positions, insertions, args.t1, args.t2)

    #if args.t is None:
    #    selected, orientation = check_positions_difference(positions, insertions)
    #else:
    #    selected, orientation = check_positions_difference(positions, insertions, args.t)


    if args.tgen is None and args.tdif is None:
        selected, orientation = check_positions_combination(positions, insertions)
    elif args.tgen is None:
        selected, orientation = check_positions_combination(positions, insertions, threshold_difference=args.tdif)
    elif args.tdif is None:
        selected, orientation = check_positions_combination(positions, insertions, threshold_genome=args.tgen)
    else:
        selected, orientation = check_positions_combination(positions, insertions, threshold_genome=args.tgen, threshold_difference=args.tdif)
    
    coverage = obtain_coverage_distribution(args.bam, insertions, selected)
    scaled = scale_distributions(coverage)
    distances = calculate_DHD(scaled)

    print_output_distance(insertions, selected, distances, args.oFilt)
    calculate_inssite(coverage, insertions, selected, orientation, args.oSite)
