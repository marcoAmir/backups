#!/usr/bin/env python
#
# 4/27/2016
# This code will intersect all the deletion regions from the species with some threshold of accuracy (by bp)

import argparse
import re
import sys
import os
import subprocess

def readArgs():
        print ""
        parser = argparse.ArgumentParser('Intersect all large deletions from different species, '
                'return a clustered deletions as bed formatted regions')
        parser.add_argument('-e', '--epsilonOverlapAccuracy', dest='epsilonOverlapAccuracy',
            help='maximum allowed number of nonoverlapping base-pairs',
            type=int, default=0)
        if len(sys.argv) > 3:
                parser.print_help()
                sys.exit(1)
        args = vars(parser.parse_args())
        return args

def run(epsilonOverlapAccuracy):
	#assemblies = ['ailMel1','bosTau7','calJac3','camFer1','canFam3','capHir1','cavPor3','cerSim1','chiLan1'
	#	,'dasNov3','dipOrd1','echTel2','eleEdw1','eptFus1','equCab2','eriEur2','felCat5'
	#	,'gorGor3','hetGla2','jacJac1','lepWed1','loxAfr3','macFas5','mesAur1','micOch1','mm10'
	#	,'musFur1','myoDav1','myoLuc2','nomLeu3','ochPri3','octDeg1','odoRosDiv1','orcOrc1'
	#	,'oryAfe1','oryCun2','otoGar3','oviAri1','panHod1','panTro4','papHam1','ponAbe2','pteAle1','pteVam1'
	#	,'rheMac3','rn5','saiBol1','sorAra2','speTri2','susScr3','tarSyr1','triMan1','tupBel1','tupChi1','turTru2','vicPac2']
	assemblies = ['bosTau7','calJac3','capHir1']
	print assemblies
	# Step 1: add all the identified deletion to a hash map
	for assembly in assemblies:
		i = len(deletionsHash)
		print i
		appendDeletions(assembly, i)
	# Step 2: add deletions from species to clusters if they meet the threshold criteria
	for assembly in assemblies:
		matchDeletionsToClusters(assembly, epsilonOverlapAccuracy)
	# Step 3: canonicalize the clusters and print output (for clusters with > 1 species deletion)
	canonicalizeClusters(epsilonOverlapAccuracy)

def appendDeletions(assembly, i):
	k = i + 1
	bedFile = assembly + '.deletions.bed'
	with open(bedFile) as f:
		for line in f:
			deletionsHash['delClust' + str(k)] = ['\t'.join(['\t'.join(line.split('\t')[:3]), assembly])]
			k += 1

def matchDeletionsToClusters(assembly, epsilonOverlapAccuracy):
	print '\n...match ' + assembly + ' deletions to large del-clusters'
	c1 = 0
	c2 = 0
	bedFile = assembly + '.deletions.bed'
	with open(bedFile) as f:
		for line in f:
			c1 += 1
			print c1
			chrom     = line.split('\t')[0]
			delStarts = int(line.split('\t')[1])
			delEnds   = int(line.split('\t')[2])
			for key, value in deletionsHash.iteritems():
				if (re.search(chrom, value[0]) and not re.search(assembly, value[0])):
					if delMatchCluster(chrom, delStarts, delEnds, value, epsilonOverlapAccuracy):
						deletionsHash[key].append('\t'.join(['\t'.join(line.split('\t')[:3]), assembly]))
						c2 += 1
	print "\t\t" + str(c2) + " deletions out of " + str(c1) + " " + assembly + " deletions were matched to large del-clusters" 

def delMatchCluster(chrom, x0, x1, value, epsilonOverlapAccuracy):
	matchCluster = 1
	for v in value:
		matchCluster = matchCluster*int(chrom==v.split('\t')[0])
		matchCluster = matchCluster*int((abs(int(v.split('\t')[1])-x0) + abs(int(v.split('\t')[2])-x1)) <= epsilonOverlapAccuracy)
	return matchCluster
		
def canonicalizeClusters(epsilonOverlapAccuracy):
	# in each cluster with more then a single deletion I'll output the region and each species in a seperate line
	# the region will span from the min of the starting point from all the deletions to the maximum of the end points
	outBed = open('overlappingDeletions.offThres_' + str(epsilonOverlapAccuracy) + 'bp.bed', 'w')
	for key, value in tmpHash.iteritems():
		if len(value)>1:
			assembliesInCluster = []
			startC = 100000000000
			endC   = -100000000000
			chrom  = 'chr1'
			for v in value:
				chrom = v.split('\t')[0]
				assembliesInCluster.append(v.split('\t')[3])
				if (int(v.split('\t')[1]) < startC):
					startC = int(v.split('\t')[1])
				if (int(v.split('\t')[2]) > endC):
					endC = int(v.split('\t')[2])
			for assembly in assembliesInCluster:
				outBed.write(chrom + '\t' + str(startC) + '\t' + str(endC) + '\t' + assembly + '\n')
	outBed.close()
			

if __name__ == "__main__":
        args = readArgs()
	deletionsHash = {}
        run(**args)

