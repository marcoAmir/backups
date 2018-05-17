#!/usr/bin/env python
#
# 3/1/2016
# This code is the module for overlapDelInProtCodRegions.py with all the relevant functions to read deletions from chain files, intersect them with sequencing gap regions etc.

import re
import sys
import os
import gzip
import subprocess

def readGapsInChainFile(assembly, subsetChains, refAssembly="hg38"):
	print("\n\t...extracting chain gaps for target assembly " + assembly + "\n")
	chainMapDir = "/cluster/u/amirma/geneLoss/" + refAssembly + "/mappings/"
	# This subroutine reads gaps in chain files (gaps of the query species seen in the reference browser)
	# and generates a BED file of those gaps with indication for whether they are single or double sided as well as chain id
	minDelSize = 1	 	# I'm interested in any non-zero deletion
	gapPad = 5	        # I'll pad the gap with gapPad in each direction to overlap with the gap track
	#chainFile='/cluster/gbdb/' + refAssembly + '/liftOver/' + refAssembly + 'To' + assembly[0].upper() + assembly[1:] + '.over.chain.gz'
	chainFile='/cluster/u/yatisht/aligner/comparative/' + refAssembly + '/chains/' + refAssembly + '.' + assembly + '.all.chain'
	#if (assembly=='orcOrc1'):
	#	chainFile='/cluster/u/amirma/data/chainHelpers/hg38/' + assembly + '/hg38.' + assembly + '.all.chain.gz'
	#chainFile = 'hg19To' + assembly[0].upper() + assembly[1:] + '.over.chain.sample.gz'
	outBed = open('/cluster/u/amirma/geneLoss/hg38/browser_visuals/largeDel/' + assembly + '.gapsInChains.bed', 'w')
	refChr = 'chr1' # just to pass the first control
	queryGapTrack = getQueryGapTrack(assembly)			# A hash table
	chainIDs = GetSubsetChains(assembly, chainMapDir, chainFile)	# an array with chain IDs (human to query to analyze) 
	#with gzip.open(chainFile) as f:
	with open(chainFile) as f:
		for line in f:
			if ((len(line.split('\t'))<=1 and len(line.split(' '))<=1) or line[0]=='#'):
				continue
			if re.search('chain', line):	# new chain
				refChr  = line.split(' ')[2]		# chr in reference genome
				chainID = line.split(' ')[12].rstrip()	# the id of the chain
				if (not conventionalChromosomeHG(refChr)) or (chainID not in chainIDs):
					while line != "\n":
						line = next(f)
					continue
				qChr    = line.split(' ')[7]		# chr in query
				qStrand = line.split(' ')[9]		# = or - strand in query
				refPos  = int(line.split(' ')[5])	# the start pos of the chain - on reference
				qPos    = int(line.split(' ')[10])	# the start pos of the chain - on query
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
							outBed.write(refChr + '\t' + str(refPos) + '\t' + str(RefEndDel) + '\t' + chainID + '_DS_gap' + str(gapInd) + '\n')
						else:
							outBed.write(refChr + '\t' + str(refPos) + '\t' + str(RefEndDel) + '\t' + chainID + '_SS_gap' + str(gapInd) + '\n')
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

def GetSubsetChains(assembly, chainMapDir, chainFile, subsetChains=True):
	chainIDs = []
	if (not subsetChains):
		with gzip.open(chainFile) as f:
			for line in f:
				if re.search('chain', line):
					chainIDs.append(line.split(' ')[12].rstrip())
	else:
		with open(chainMapDir + assembly + '.chain_ids') as f:
			for line in f:
				chainIDs.append(line.split(' ')[1])
	return list(set(chainIDs))
	

def conventionalChromosomeHG(refChr):
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
			return False
	return True

def intersectDelWithCodingExons(assembly):
	print("\n\t...intersecting chain gaps with coding exons\n")
	os.system("overlapSelect ./data/testedTranscripts.ExonsCodingBases.bed ./largeDel/" + assembly + ".gapsInChains.bed stdout | sort -t$'\\t' -k4,4 > a.bed" )	
	os.system("overlapSelect ./data/testedTranscripts.ExonsCodingBases.bed ./largeDel/" + assembly + ".gapsInChains.bed -idOutput stdout | sort -u > b.txt")
	os.system("join -t$'\\t' -1 4 -2 1 -o '2.2 1.1 1.2 1.3 1.4' a.bed b.txt > ./largeDel/" + assembly + ".overlappingGaps")
	os.system("rm -rf a.bed b.txt")
	#subprocess.call(overlapSelect, shell=True)
#	if not os.path.isfile('refTranscripts.bed'):
#		transcripts = [line.strip() for line in open(refTrans, 'r')]
#		transcriptBED = open('refTranscripts.bed', 'w')
#		with open (dataDir + '/Ensembl/transcriptExonsCodingBases.bed') as f:
#			for line in f:
#				if line.split('\t')[3] in transcripts:
#					transcriptBED.write(line)
#		transcriptBED.close()
#	os.system("overlapSelect refTranscripts.bed " + assembly + ".deletions.bed " + assembly + "_tmp")
#	os.system("mv " + assembly + "_tmp " + assembly + ".deletions.bed")

			
