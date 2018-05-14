#!/usr/bin/env python
#
# 4/1/2016
# This code will iterate through intact transcript and will identify large deletions in protein-coding genes that are common to two or more adjacent species.
# the generated BED files will assign the following annotation for every entry in the 4th column:
#	assemnly_gapIndex_chainID_gapType	# e.g., mm10_1_1668438_SS the gap index means that I number the gaps along the chain (as appear in hg19).
						# gapType: SS means single sided gap; DS_2 means double sided gap and show the gap size in the query species 

import argparse
import os
import re
import socket
import sys

import readLargeDeletions as RDEL

def readArgs():
	print ""
	parser = argparse.ArgumentParser('Read large deletions in query species, interesct with protein-coding exons, '
		'and identify overlapping deletions in 2 or more adjacent taxa')
	parser.add_argument('assembly', metavar='assembly',
		help='The accurate name+number of the query species assembly')
	parser.add_argument('refTrans', metavar='refTranscript.txt',
		help='Either a file with a list of Ensembl transcript ID or a single ENST identifier')
	parser.add_argument('dataDir', metavar='$GENELOSS/data/',
		help='Path to directory with gene losses screen data: /cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data')
	parser.add_argument('-s', '--subsetChains', dest='subsetChains',
            help='whether or not (1/0) to read only from chains that uniquely map transcripts',
            type=int, default=1)
	parser.add_argument('-g', '--intersectGenes', dest='intersectGenes',
            help='whether or not (1/0) to leave only gaps that intersect coding exons',
            type=int, default=1)
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	return args

def run(assembly, refTrans, dataDir, subsetChains, intersectGenes):
	RDEL.readGapsInChainFile(assembly, subsetChains, dataDir)
	if (intersectGenes!=0):
		RDEL.intersectDelWithCodingExons(assembly, refTrans, dataDir)

if __name__ == "__main__":
        args = readArgs()
	run(**args)
