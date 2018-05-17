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

import readDeletionsFromChains as RDEL

transcript_coding_exons = "/cluster/u/amirma/geneLoss/hg38/browser_visuals/data/testedTranscripts.ExonsCodingBases.bed"

def readArgs():
	print ""
	parser = argparse.ArgumentParser('Read large deletions in query species, and intersects with protein-coding exons')
	parser.add_argument('assembly', metavar='assembly',
		help='The accurate name+number of the query species assembly')
	parser.add_argument('-s', '--subsetChains', dest='subsetChains',
            help='whether or not to read only from chains that uniquely map transcripts',
            type=bool, default=True)
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	return args

def run(assembly, subsetChains):
	#RDEL.readGapsInChainFile(assembly, subsetChains)
	RDEL.intersectDelWithCodingExons(assembly)


if __name__ == "__main__":
        args = readArgs()
	run(**args)
