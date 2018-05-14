#!/usr/bin/env python
#
# 4/1/2016
# This code is the module for overlapDelInProtCodRegions.py with all the relevant functions to read deletions from chain files, intersect them with sequencing gap regions etc.

import re
import sys
import os
import gzip
import subprocess

def readGapsInChainFile(assembly, subsetChains, dataDir):
	# This subroutine reads gaps in chain files (gaps of the query species seen in the hg19 browser)
	# and generates a BED file of those gaps with indication for whether they are single or double sided as well as chain id
	minDelSize = 20 	# I'll default the deletion size of interest to this value
	gapPad = 5	        # I'll pad the gap with gapPad in each direction to overlap with the gap track
	chainFile='/cluster/gbdb/hg19/liftOver/hg19To' + assembly[0].upper() + assembly[1:] + '.over.chain.gz'
	if (assembly=='orcOrc1'):
		chainFile='/cluster/u/amirma/data/chainHelpers/hg19/' + assembly + '/hg19.' + assembly + '.all.chain.gz'
	#chainFile = 'hg19To' + assembly[0].upper() + assembly[1:] + '.over.chain.sample.gz'
	outBed = open(assembly + '.deletions.bed', 'w')
	refChr = 'chr1' # just to pass the first control
	queryGapTrack = getQueryGapTrack(assembly)			    # A hash table
	chainIDs = GetSubsetChains(assembly, subsetChains, dataDir, chainFile) # an array with chain IDs (human to query to analyze) 
	with gzip.open(chainFile) as f:
		for line in f:
			if ((len(line.split('\t'))<=1 and len(line.split(' '))<=1) or line[0]=='#'):
				continue
			if re.search('chain', line):	# new chain
				refChr  = line.split(' ')[2]		# chr in hg19
				chainID = line.split(' ')[12].rstrip()	# the id of the chain
				if (not conventionalChromosomeHG19(refChr)) or (chainID not in chainIDs):
					while line != "\n":
						line = next(f)
					continue
				#else:
				#	print line
				qChr    = line.split(' ')[7]		# chr in query
				qStrand = line.split(' ')[9]		# = or - strand in query
				refPos  = int(line.split(' ')[5])	# the start pos of the chain
				qPos    = int(line.split(' ')[10])	# the start pos of the chain
				gapInd  = 0				# index of the gap
				if (qStrand=="+"):			# I'll coerce the query strand to an actual direction for subsequent operation on genomic regions
					qStrandM = 1
				if (qStrand=="-"):
					qStrandM = -1
					qPos = int(line.split(' ')[8]) - qPos
			else:
				refPos = refPos + int(line.split('\t')[0])	# Propagate the ref position to the end of current alignment block
				qPos = qPos + qStrandM*int(line.split('\t')[0])
				refDelSize = int(line.split('\t')[1])
				qDelSize   = int(line.split('\t')[2])
				RefEndDel  = refPos + refDelSize
				if (refDelSize>=minDelSize):
					checkGapInterval = [qPos - qStrandM*gapPad, qPos + qStrandM*(gapPad + qDelSize)]
					checkGapInterval.sort()
					if (qChr not in queryGapTrack.keys()) or (notOnSeqeuncingGap(checkGapInterval, queryGapTrack[qChr])):	# Check if in seqeuncing gap
						gapInd += 1
						if (qDelSize>0):
							outBed.write(refChr + '\t' + str(refPos) + '\t' + str(RefEndDel) + '\t' + assembly + '_' + str(gapInd) + '_' + chainID + '_DS_' + str(qDelSize) + '\n')
						else:
							outBed.write(refChr + '\t' + str(refPos) + '\t' + str(RefEndDel) + '\t' + assembly + '_' + str(gapInd) + '_' + chainID + '_SS' + '\n')
				refPos = RefEndDel
				qPos = qPos + qStrandM*qDelSize
	outBed.close()

def getQueryGapTrack(assembly):
	gapFile = '/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data/UCSC/GapTracks/' + assembly + '.gap_filtered.bed'
	gapHash = {}
	with open(gapFile) as g:
		for line in g:
			chrom    = line.split('\t')[0]
			gapStart = line.split('\t')[1]
			gapEnd   = line.split('\t')[2]
			if chrom in gapHash.keys():
				gapHash[chrom].append(gapStart + '\t' + gapEnd)
			else:
				gapHash[chrom] = [gapStart + '\t' + gapEnd]
	return gapHash

def GetSubsetChains(assembly, subsetChains, dataDir, chainFile):
	chainIDs = []
	if (subsetChains==0):
		with gzip.open(chainFile) as f:
			for line in f:
				if re.search('chain', line):
					chainIDs.append(line.split(' ')[12].rstrip())
	else:
		with open(dataDir + '/hg19.calltable.unique.placentals') as f:
			for line in f:
				if (line.split('\t')[0]==assembly):
					chainIDs.append(line.split('\t')[2].split('_')[0])
	return list(set(chainIDs))
	

def conventionalChromosomeHG19(refChr):
	a = (re.search('chr\d+\Z', refChr) and int(re.search('chr(\d+)\Z', refChr).group(1))<=22 and int(re.search('chr(\d+)\Z', refChr).group(1))>=1)
	b = (refChr=='chrM') or (refChr=='chrX') or (refChr=='chrY')
	if (a or b):
		return 1
	else:
		return 0
		
	
def notOnSeqeuncingGap(chainGap, ChromQueryGapTrack):
	# I check if a gap in chain contains or has any overlap with a sequencing gap
	x0 = chainGap[0]	
	y0 = chainGap[1]
	l0 = y0 - x0	# total length of gap in chain
	for seqGap in ChromQueryGapTrack:
		x1 = int(seqGap.split('\t')[0])
		y1 = int(seqGap.split('\t')[1])
		l1 = y1 - x1	# total length of current sequencing gap
		l_observed = max(y0, y1) - min(x0, x1)
		if (l_observed < (l0 + l1)):
			return 0
	return 1

def intersectDelWithCodingExons(assembly, refTrans, dataDir):
	if not os.path.isfile('refTranscripts.bed'):
		transcripts = [line.strip() for line in open(refTrans, 'r')]
		transcriptBED = open('refTranscripts.bed', 'w')
		with open (dataDir + '/Ensembl/transcriptExonsCodingBases.bed') as f:
			for line in f:
				if line.split('\t')[3] in transcripts:
					transcriptBED.write(line)
		transcriptBED.close()
	os.system("overlapSelect refTranscripts.bed " + assembly + ".deletions.bed " + assembly + "_tmp")
	os.system("mv " + assembly + "_tmp " + assembly + ".deletions.bed")

			
