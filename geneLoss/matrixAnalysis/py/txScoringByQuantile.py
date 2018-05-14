#!/usr/bin/env python
#
# 9/2/2015
# This code reads an input of pairwise alignments attributes (Transcript-species alignment with counts of frameshifts, stop codons, deletions), and reassign a LoF call (1) for a transcript in a species, based on quantiles

import argparse
import os
import re
import socket
import sys

import numpy as np

def readArgs():
        print ""
        parser = argparse.ArgumentParser('Reassignments of LoF calls from quantiles')
        parser.add_argument('callTable', metavar='$GENELOSS/data/hg19.calltable.unique.placentals',
                help='Path to a call table file, e.g., : /cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data/hg19.calltable.unique.placentals')
	parser.add_argument('dataDir', metavar='$GENELOSS/data/',
		help='Path to directory with gene losses screen data: /cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data')
	parser.add_argument('output', metavar='quantile.out',
		help='A name for the output file')
	parser.add_argument('-quantThrs', '--quantThrs', dest='quantThrs',
		help='Specify the quantile threshold for calling transcript LoF',
            	type=float, default=90.0)
        if len(sys.argv) < 4:
                parser.print_help()
                sys.exit(1)
        args = vars(parser.parse_args())
        return args

def run(callTable, dataDir, output, quantThrs):
	speciesQuantile = readSpeciesQuantiles(dataDir + "speciesList.txt", dataDir + "quantiles/")
	speciesTxKaKs = getKaKs(dataDir + "speciesList.txt", dataDir + "KaKs/speciesData/")
	o = open(output, 'w')
	f = open(callTable)
	o.write(f.readline())		# read the header to the output
	#c = 1
	for line in f:
		#print c
		#c += 1
		KaKs = speciesTxKaKs[line.split()[0] + "_" + line.split()[1]]
		mq = meanQuantile(speciesQuantile[line.split()[0]], line, KaKs)
		o.write("\t".join(line.split()[0:23]))
		if(mq > quantThrs):
			o.write("\t%d\n" %(1))
		else:
			o.write("\t%d\n" %(0))
	f.close()
	o.close()

def readSpeciesQuantiles(speciesFile, quantileDir):
	assemblies = {}
	f = open(speciesFile)
	for line in f:
		if (line.split()[2]=="primate" or line.split()[2]=="rodent" or line.split()[2]=="placental"):
			spQuantile = np.genfromtxt(quantileDir + line.split()[3] + ".canonical.Quantiles", skip_header=1)
			assemblies[line.split()[3]] = spQuantile
	f.close()
	return assemblies

def getKaKs(speciesFile, KaKsDir):
	KaKs = {}	# A hash map that will hold KaKs value for a given transcript-species pairwise, where keys are in the form: mm10_ENST00000002501
	f = open(speciesFile)
	for line1 in f:
		if (line1.split()[2]=="primate" or line1.split()[2]=="rodent" or line1.split()[2]=="placental"):
			readAssembly = line1.split()[3]
			g = open(KaKsDir + "transcriptNonSynRate." + readAssembly + ".tsv")
			next(g)
			for line2 in g:
				key = readAssembly + "_" + line2.split()[1]
				KaKs[key] = line2.split()[12]
			g.close()
	f.close()
	return KaKs
			

def meanQuantile(spQuantile, TxPwAttr, KaKs):
	fracDelExons = float(TxPwAttr.split()[4])/float(TxPwAttr.split()[14])
	fracDelBases = float(TxPwAttr.split()[9])/float(TxPwAttr.split()[13])
	quantileStop         = quantileFromColumn(float(TxPwAttr.split()[3]), spQuantile[:,[0,1]])
	quantileFrameshifts  = quantileFromColumn(float(TxPwAttr.split()[10]), spQuantile[:,[0,2]])
	quantileDelExons     = quantileFromColumn(float(TxPwAttr.split()[4]), spQuantile[:,[0,3]])
	quantileDelBases     = quantileFromColumn(float(TxPwAttr.split()[9]), spQuantile[:,[0,5]])
	quantileFracDelExons = quantileFromColumn(float(fracDelExons), spQuantile[:,[0,4]])
	quantileFracDelBases = quantileFromColumn(float(fracDelBases), spQuantile[:,[0,6]])
	if(type(KaKs) == float):
		quantileKaKs = quantileFromColumn(KaKs, spQuantile[:,[0,7]])
		return np.mean([quantileStop, quantileFrameshifts, quantileDelExons, quantileDelBases, quantileFracDelExons, quantileFracDelBases, quantileKaKs])
	else:
		return np.mean([quantileStop, quantileFrameshifts, quantileDelExons, quantileDelBases, quantileFracDelExons, quantileFracDelBases])
	
def quantileFromColumn(x, u):
	for i in range(0,len(u)):
		if (x == u[i,1]):
			return u[i,0]
		if (x < u[i,1]):
			return u[i-1,0]
	return 100.0
	

if __name__ == "__main__":
        args = readArgs()        
	run(**args)
