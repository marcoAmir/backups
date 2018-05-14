#!/usr/bin/env python
#
# 9/2/2015
# This code reads an input species/assembly and a given reference transcript(s) and returns for each transcript on which exons does he see frameshifts

import argparse
import os
import re
import socket
import sys

def readArgs():
	print ""
	parser = argparse.ArgumentParser('Reads non 0 modulu 3 exon indels in a query species'
		', for a given assembly and transcript ID')
	parser.add_argument('assembly', metavar='query',
		help='Provide a species or assembly (currently supports the latter)')
	parser.add_argument('refTrans', metavar='refTranscript.txt',
		help='Either a file with a list of Ensembl transcript ID or a single ENST identifier')
	parser.add_argument('dataDir', metavar='$GENELOSS/data/',
		help='Path to directory with gene losses screen data: /cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data')
	if len(sys.argv) < 3:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	return args

def run(assembly, refTrans, dataDir):
	print "\nNote that with the current version, reference exons overlapping any sequencing gaps are ignored in the count\n"
	inFile = dataDir+"output/hg19."+assembly+".codingMutations.out"
	outFile = assembly+".hg19transcripts.statusByExons.txt"
	transcriptHash = getTranscriptsBed(dataDir)
	print "\n\t...Parsing "+inFile+"\n" 
	enstQueryMap = {}
	with open(inFile) as f:
		for line in f:
			enstQueryMap[line.split()[1]] = line.split()[8]
	f.close()
	if re.match( '^ENST\d+$', refTrans):
		writeFrameshifts(assembly, enstQueryMap, transcriptHash, refTrans, outFile, 0)
	else:
		f = open(refTrans)
		out = open(outFile, 'w')
		out.write("assembly	enstID	exon	frameshifts\n")
		out.close()
		for line in f:
			currTrans = line.rstrip()
			writeFrameshifts(assembly, enstQueryMap, transcriptHash, currTrans, outFile, 1)
		f.close() 

def getTranscriptsBed(dataDir):
	bed = open(dataDir+"/Ensembl/transcriptExonsCodingBases.bed")
	hashMap = {}
	for line in bed:
		hashMap[line.split()[3]] = [line.split()[0], line.split()[1], line.split()[2], line.split()[10], line.split()[11], line.split()[5]]
	bed.close()
	return hashMap
		
def writeFrameshifts(assembly, enstQueryMap, transcriptHash, refTrans, outFile, multiTrans):
	# The last argument (getting 1/0) will indicate whether the input is a single transcript or multiple transcript IDs in an input file
	if multiTrans:
		out = open(outFile, 'a')	
	indelRegions = open('indels.'+assembly+'.loci', 'a')
	exonDel = open('exonDeletions.'+assembly+'.loci', 'a')
	totalFramshifts = 0
	print refTrans
	for i in range(len(enstQueryMap[refTrans].split("|"))):
		if ('?' in enstQueryMap[refTrans].split("|")[i]) and multiTrans:
			out.write("%s\t%s\t%d\tNA\n" %(assembly, refTrans, i+1))
			print "query sequence gap in exon ",i+1
		#if ('-' in enstQueryMap[refTrans].split("|")[i]):
		else:
			frameshifts = checkIndel(enstQueryMap[refTrans].split("|")[i])
			if (frameshifts>0):
				print frameshifts," frameshift identified in exon ",i+1
				totalFramshifts += frameshifts
				bedCoords = getFrameshiftsLoci(i+1, enstQueryMap[refTrans].split("|")[i], transcriptHash[refTrans])
				for b in bedCoords:
					indelRegions.write(b+"\t"+refTrans+"\texon"+str(i+1)+"\n")
			if (frameshifts == -1):
				print "likely a complete deletion in exon ",i+1
				bedCoords = getFrameshiftsLoci(i+1, enstQueryMap[refTrans].split("|")[i], transcriptHash[refTrans])
				for b in bedCoords:
					exonDel.write(b+"\t"+refTrans+"\texon"+str(i+1)+"\n")
			if (frameshifts == 0):
				print frameshifts," indels identified in exon ",i+1
			if multiTrans:
				out.write("%s\t%s\t%d\t%d\n" %(assembly, refTrans, i+1, frameshifts))
			
	print "====>",totalFramshifts," TOTAL FRAMESHIFTS IDENTIFIED IN",refTrans,"IN",assembly,"\n"
	if multiTrans:
		out.close()
	indelRegions.close()
	exonDel.close()

		
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
    
if __name__ == "__main__":	
	args = readArgs()
	run(**args)
	
