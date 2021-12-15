#!/usr/bin/env python3

def UmiToSeq(inputFile, outputFile):
    fasta = SeqIO.parse(inputFile, "fasta")
    out = open(outputFile, 'w')
    params = {}

    for read in fasta:
        paramsList = read.id.split(";")
        params[paramsList[0]] = {}
        for param in paramsList[1:]:
            params[paramsList[0]][param.split("=")[0]] = param.split("=")[1]

    for k in params.keys():
        out.write(">{}\n{}\n".format(k, params[k]["seq"]))

    out.close()




if __name__ == "__main__":
    from Bio import SeqIO
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="Fasta file containing UMIs", action="store", required=True)
    parser.add_argument("-o", help="Output fasta file to contain sequences", action="store", required=True)
    args = parser.parse_args()

    UmiToSeq(args.i, args.o)
    
