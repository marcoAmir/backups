#!/usr/bin/env python
#
# 10/20/2015
# This code will iterate through intact transcript and will identify genomic events that are common to two adjacent species.

import argparse
import os
import re
import socket
import sys

def readArgs():
	print ""
	parser = argparse.ArgumentParser('Write a summary of all the events (fs, stop, deletion) per transcript'
		', such that likely ancestral events (shared between two adjacent taxa) are easily identified')
	parser.add_argument('refTrans', metavar='refTranscript.txt',
		help='Either a file with a list of Ensembl transcript ID or a single ENST identifier')
	parser.add_argument('dataDir', metavar='$GENELOSS/data/',
		help='Path to directory with gene losses screen data: /cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data')
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	return args

def run(refTrans, dataDir):
	assembliesTransPairs, assembliesTransAminoAcids = getAssembliesTrans(dataDir, refTrans)
	transcriptHash = getTranscriptsBed(dataDir)
	frameshiftHash = {}
	delHash = {}
	insHash = {}
	stopCodonHash = {}
	out = open('t', 'w')	#toDel
	for currTxAssembly in assembliesTransPairs:
		print currTxAssembly
		pairwise = assembliesTransPairs[currTxAssembly]			# the status of transcript in pairwise alignment by coding exons
		protSequence = assembliesTransAminoAcids[currTxAssembly]	# the amino-acid sequence of the transcript in the query
		if re.search('\*', protSequence):
			stopCodonHash = checkStopCodons(currTxAssembly, protSequence)	# if any stop codons were found in a particular species, returns array of loci on hg19 
		for n in range(0, len(pairwise.split('|'))):
			if ('?' in pairwise.split('|')[n]):
				continue
			frameshifts = checkIndel(pairwise.split('|')[n])
			if (frameshifts>0):
				bedCoords = getFrameshiftsLoci(n+1, pairwise.split('|')[n], transcriptHash[currTxAssembly.split("_")[1]])
				bedCoordDel, bedCoordIns = getFrameshiftsLoci2(n+1, pairwise.split('|')[n], transcriptHash[currTxAssembly.split("_")[1]])
				i = 1
				for loci in bedCoords:
					out.write(currTxAssembly + "_" + str(i) + "\texon" + str(n+1) + "\t" + loci + "\n")
					frameshiftHash[currTxAssembly + "_exon" + str(n+1) + "_" + str(i)] = "exon" + str(n+1) + "\t" + loci
					i += 1
				i = 1
				for loci in bedCoordDel:
					delHash[currTxAssembly + "_exon" + str(n+1) + "_" + str(i)] = "exon" + str(n+1) + "\t" + loci
					i += 1
				i = 1
				for loci in bedCoordIns:
					insHash[currTxAssembly + "_exon" + str(n+1) + "_" + str(i)] = "exon" + str(n+1) + "\t" + loci
					i += 1
	out.close() #toDel	
	printMutations(frameshiftHash, 'frameshift', 'frameshifts.txt')
	# NOW , FIND ALL THE KEYS WITH SIMILAR VALUES OF FRAMESHIFTS - good for QA purpose, maybe I'll remove these parts later on
	readSharedFrameshifts(frameshiftHash, 'frameshifts.common')
	readSharedFrameshifts(insHash, 'insertions.common')
	readSharedFrameshifts(delHash, 'deletions.common')
	# AND GENERATE SIMILAR MAPPPINGS FOR EXON DELETIONS AND STOP-CODONS
	readSharedStopcodons(stopCodonHash, 'stopCodons.common')
	

def getAssembliesTrans(dataDir, refTrans):
	transcripts = []
	f1 = open(refTrans)
	for line in f1:
		transcripts.append(line.strip())
	f1.close()
	assembliesTransPairs = {}
	assembliesTransAminoAcids = {}
	f2 = open(dataDir + "/hg19.calltable.unique.placentals")
	f2.next()
	assembly = 'NA'
	for line in f2:
		if (line.split()[0] != assembly):
			assembly = line.split()[0]
			print ('\n\t...Reading ' + assembly)
			TxStatus = getTransStatInSpecies(dataDir, assembly)
		#if (line.split()[1] in transcripts):
		key = line.split()[0] + "_" + line.split()[1]
		assembliesTransPairs[key]      = TxStatus[line.split()[1]][0]	# indels data by exons
		assembliesTransAminoAcids[key] = TxStatus[line.split()[1]][1]	# amino acid sequence in query species
	f2.close()
	return assembliesTransPairs, assembliesTransAminoAcids

def getTranscriptsBed(dataDir):
	bed = open(dataDir + "/Ensembl/transcriptExonsCodingBases.bed")
	hashMap = {}
	for line in bed:
		hashMap[line.split()[3]] = [line.split()[0], line.split()[1], line.split()[2], line.split()[10], line.split()[11], line.split()[5]]
	bed.close()
	return hashMap

def getTransStatInSpecies(dataDir, assembly):
	TxStatus = {}
	f = open(dataDir + "/output/hg19." + assembly + ".codingMutations.out")
	for line in f:
		TxStatus[line.split()[1]] = line.split()[8:10]
	f.close()
	return TxStatus

def checkStopCodons(currTxAssembly, protSequence):
	s_out = open('stopCodons.QA', 'a')
	protSequenceProc = ''.join([i for i in protSequence if not i.isdigit()])	# this just removes counts of bp substitutions
	transcript = transcriptHash[currTxAssembly.split("_")[1]]
	blocks = transcript[3].split(",")
	alignment = assembliesTransPairs[currTxAssembly].split("|")
	aminoAcidsByExons = []
	# rearrange the exons by mutiples of 3 (an excess of 1bp moves to the next blocks; an excess of 2 bases - i.e., one is missing, we take 1 bp from the next block)
	blocksMod3 = transcript[3].split(",")
	for i in range(0,len(blocksMod3)-1):
		if(int(blocksMod3[i])%3 == 1):
			blocksMod3[i]   = str(int(blocksMod3[i])-1)
			blocksMod3[i+1] = str(int(blocksMod3[i+1])+1)
		if(int(blocksMod3[i])%3 == 2):
			blocksMod3[i]   = str(int(blocksMod3[i])+1)
			blocksMod3[i+1] = str(int(blocksMod3[i+1])-1)
	# Step 1: First I will split the amino acid sequence by exons:
	for i in range(0, len(alignment)):
		if ((not ("-" in alignment[i])) or ("0-" in alignment[i])) and (not re.search('[1-9]+\-', alignment[i])):	# no insertion, take the number of amino acids as the length of block/3
			aminoAcidsByExons.append(protSequenceProc[0:int(blocksMod3[i])/3])
			protSequenceProc = protSequenceProc[int(blocksMod3[i])/3:]
		else:
			l = 0	# block size in the query
			for alignmentBlock in alignment[i].split(","):
				l += int(alignmentBlock.split("-")[0].replace("?",""))
			if(l%3 == 0):
				aminoAcidsByExons.append(protSequenceProc[0:l/3])
				protSequenceProc = protSequenceProc[l/3:]
			else:	# if I find 1 or 2 excess bp I'll add a '?' 
				aminoAcidsByExons.append(protSequenceProc[0:l/3] + "?")
				protSequenceProc = protSequenceProc[l/3:]
		#s_out.write(currTxAssembly + "\t" + blocks[i] + "\t" + alignment[i] + "\t" + aminoAcidsByExons[len(aminoAcidsByExons)-1] + "\t" + str(float(len(aminoAcidsByExons[len(aminoAcidsByExons)-1]))*3) + "\n")
		print(currTxAssembly + "\t" + blocks[i] + "\t" + blocksMod3[i] + "\t" + alignment[i] + "\t" + aminoAcidsByExons[len(aminoAcidsByExons)-1] + "\t" + str(float(len(aminoAcidsByExons[len(aminoAcidsByExons)-1]))*3) + "\n")
	#s_out.write(currTxAssembly + "\t" + "leftover:\t" + protSequenceProc + "\t\t" + "length: " + str(len(protSequenceProc)) + "\n\n")
	
	# Step 2: Now, I will iterate exon by exon, and if I found a stop codon I'll call a function that extracts its hg19 coordinates:
	locis = []
	for i in range(0, len(alignment)):
		if re.search('\*', aminoAcidsByExons[i]):
			flanks=[]	# for QAing purposes I'll register the flanking 2 amino-acids to the stop-codon
			for p in [pos for pos, char in enumerate(aminoAcidsByExons[i]) if char == '*']:
				if (p-3 < 0):
					flanks.append(aminoAcidsByExons[i][:p+4])
				else:
					flanks.append(aminoAcidsByExons[i][p-3:p+4])
			locis = getStopOnHg19(transcript, i, aminoAcidsByExons[i])	# call a function that computes the location of the stop codon on hg19 (need to write it)
			for j in range(0,len(locis)):
				s_out.write(currTxAssembly + "\t" + locis[j] + "\t" + flanks[j] + "\t" + transcript[5] +"\n")
				stopCodonHash[currTxAssembly.split("_")[0] + "|" + flanks[j]] = currTxAssembly.split("_")[1] + "\t" + locis[j]
	s_out.close()
	return stopCodonHash
	#return locis
			
def getStopOnHg19(transcript, i, aminoAcids):
	locis = []
	if (transcript[5] == "+"):	# + strand
		start = int(transcript[4].split(",")[i]) + int(transcript[1])
		s = 0
		for j in range(0, len(aminoAcids)):
			if re.match('\*', aminoAcids[j]):
				tickBed = transcript[0] + "\t" + str(start + s) + "\t" + str(start + s + 2)
				s += 3
				locis.append(tickBed)
			else:
				s += 3
	else:				# - strand
		start = int(transcript[4].split(",")[i]) + int(transcript[3].split(",")[i]) + int(transcript[1])
		s = 0
		for j in range(0, len(aminoAcids)):
			if re.match('\*', aminoAcids[j]):
				tickBed = transcript[0] + "\t" + str(start - s - 2) + "\t" + str(start - s)
				s += 3
				locis.append(tickBed)
			else:
				s += 3
	return locis

				

	
	

def checkIndel(exon):
	# This subroutine will take an exon alignment block (a1,b1,a2,b2...) and will determine how many non 0 modulu 3 indels it has
	# The function can return either:
	# (1) number of non zero modulu 3 indels (an integer > 0)
	# (2) 0 if non indel detected or if indels are 0 modulu 3
	# (3) -1 if it is likely that all or most of the exon is deleted (alignment block of length 0 at the begining and end) 
	indelPerExon = 0
	if int(exon.split(",")[0])==0 and int(exon.split(",")[len(exon.split(","))-1])==0:
		return -1
	else:
		for alignmentBlock in exon.split(","):
			if ('-' in alignmentBlock) and (exon.split(",")[0]>0) and (exon.split(",")[len(exon.split(","))-1]>0):
				indel = int(alignmentBlock.split("-")[0])-int(alignmentBlock.split("-")[1])
				if (abs(indel)%3):
					indelPerExon += 1
	return indelPerExon

def getFrameshiftsLoci(n ,exon, bed):
	#print exonSize, exonStart
	locis = []
	blocks = bed[3]
	starts = bed[4]
	strand = bed[5]
	nExons = len(blocks.split(",")) - 1
	if(strand == "+"):
		exonSize = blocks.split(",")[n-1]
		exonStart = int(bed[1]) + int(starts.split(",")[n-1])
		l = exonStart	# here I'll sum the aligning blocks	
		for alignmentBlock in exon.split(","):
			if ('-' in alignmentBlock):
				insertion = alignmentBlock.split("-")[0]
				deletion = alignmentBlock.split("-")[1]
				if (abs(int(insertion)-int(deletion))%3):
					tick = int(deletion) + int(l)
					tickBed = bed[0]+"\t"+str(l)+"\t"+str(tick)
					locis.append(tickBed)
					l = int(tick)
			else:
				l += int(alignmentBlock)
	if(strand == "-"):
		exonSize = blocks.split(",")[nExons-n]
		exonStart = int(bed[1]) + int(starts.split(",")[nExons-n]) + int(exonSize) 
		l = exonStart	# here I'll sum the aligning blocks	
		for alignmentBlock in exon.split(","):
			if ('-' in alignmentBlock):
				insertion = alignmentBlock.split("-")[0]
				deletion = alignmentBlock.split("-")[1]
				if (abs(int(insertion)-int(deletion))%3):
					tick = int(l) - int(deletion)
					tickBed = bed[0]+"\t"+str(tick)+"\t"+str(l)
					locis.append(tickBed)
					l = int(tick)
			else:
				l -= int(alignmentBlock)
	return locis

def getFrameshiftsLoci2(n ,exon, bed):
	#print exonSize, exonStart
	locisDel = []
	locisIns = []
	blocks = bed[3]
	starts = bed[4]
	strand = bed[5]
	nExons = len(blocks.split(",")) - 1
	if(strand == "+"):
		exonSize = blocks.split(",")[n-1]
		exonStart = int(bed[1]) + int(starts.split(",")[n-1])
		l = exonStart	# here I'll sum the aligning blocks	
		for alignmentBlock in exon.split(","):
			if ('-' in alignmentBlock):
				insertion = int(alignmentBlock.split("-")[0])
				deletion = int(alignmentBlock.split("-")[1])
				if (abs(insertion-deletion)%3):
					tick = deletion + int(l)
					tickBed = bed[0]+"\t"+str(l)+"\t"+str(tick)
					if (insertion>0 and deletion==0):
						locisIns.append(tickBed)
					if (insertion==0 and deletion>0):
						locisDel.append(tickBed)
					if (insertion>0 and deletion>0):
						locisIns.append(tickBed)
						locisDel.append(tickBed)
					l = int(tick)
			else:
				l += int(alignmentBlock)
	if(strand == "-"):
		exonSize = blocks.split(",")[nExons-n]
		exonStart = int(bed[1]) + int(starts.split(",")[nExons-n]) + int(exonSize) 
		l = exonStart	# here I'll sum the aligning blocks	
		for alignmentBlock in exon.split(","):
			if ('-' in alignmentBlock):
				insertion = int(alignmentBlock.split("-")[0])
				deletion = int(alignmentBlock.split("-")[1])
				if (abs(insertion-deletion)%3):
					tick = int(l) - deletion
					tickBed = bed[0]+"\t"+str(tick)+"\t"+str(l)
					if (insertion>0 and deletion==0):
						locisIns.append(tickBed)
					if (insertion==0 and deletion>0):
						locisDel.append(tickBed)
					if (insertion>0 and deletion>0):
						locisIns.append(tickBed)
						locisDel.append(tickBed)
					l = int(tick)
			else:
				l -= int(alignmentBlock)
	return (locisDel, locisIns)

def printMutations(mutationHash, mutationType, filename):
	outMut = open(filename, 'w')
	if mutationType=='frameshift':
		for k in mutationHash:
			species  = ''.join([i for i in k.split("_")[0] if not i.isdigit()]) 	# assembly
			mutation = '_'.join([k.split("_")[1], '_'.join(frameshiftHash[k].split('\t')[1:])])  # loci
			outMut.write(mutation + "\t" + species + "\t1\n")	
	outMut.close()	

def readSharedFrameshifts(frameshiftHash, filename):
	invframeshiftHash = {}
	for k, v in frameshiftHash.iteritems():
		invframeshiftHash.setdefault(k.split("_")[1] + "\t" + v, []).append(k.split("_")[0])
	outFs = open(filename, 'w')
	for bed in invframeshiftHash:
		species = list(set(invframeshiftHash[bed]))
		species.sort()
		outFs.write(bed + "\t" + ",".join(species) + "\t" + str(len(species)) + "\n")
	outFs.close()

def readSharedStopcodons(stopCodonHash, filename):
	invstopCodonHash = {}
	for k, v in stopCodonHash.iteritems():
		invstopCodonHash.setdefault(v, []).append(k)
	outStops = open(filename, 'w')
	for bed in invstopCodonHash:
		species = list(set(invstopCodonHash[bed]))
		species.sort()
		outStops.write(bed + "\t" + ",".join(species) + "\t" + str(len(species)) + "\n")
	outStops.close

if __name__ == "__main__":
        args = readArgs()        
	run(**args)
