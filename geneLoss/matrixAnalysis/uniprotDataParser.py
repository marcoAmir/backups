#!/usr/bin/env python
#
# 10/30/2015
# This code reads data from Uniprot-swissprot (reviewed gene-protein models) and parse it into a tab delimited file (1 row per entry) with the following columns:
#	col1 - uniprot ID. (I've selected the primary one)	example: P17718
#	col2 - the secondary identifier (the parenthesized)	example: 2SS_MATST
#	col3 - protein length in amino-acids
#	col4 - taxon ID						example: 9606
#	col5 - species name					example: Bos taurus (Bovine).
#	col6 - taxonomical classification			example: Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Laurasiatheria; Cetartiodactyla; Ruminantia; Pecora; Bovidae; Bovinae; Bos.
#	col7 - protein evidence level				example:  1: Evidence at protein level;

import argparse
import os
import re
import socket
import sys
import gzip

def readArgs():
        print ""
        parser = argparse.ArgumentParser('A parser for Uniprot/Swissprot database')
        parser.add_argument('uniprotData', metavar='uniprot_sprot.dat.gz',
                help='A zipped file with the uniprot data')
        if len(sys.argv) < 1:
                parser.print_help()
                sys.exit(1)
        args = vars(parser.parse_args())
        return args

def run(uniprotData):
	uniEntry = {}
	lineage = ''
	species = ''	
	ID1_empty = 1
	f = gzip.open(uniprotData, 'r')
	c = 1
	for line in f:
		print c
		c = c + 1
		currLine = line.rstrip()
		if re.search('//$', currLine):
			printEntry(uniEntry)
			uniEntry = {}
			uniEntry['EnsID'] = ''
			lineage = ''
			species = ''
			ID1_empty = 1
		else:
			if re.split('\s+', currLine)[0]=="ID":
				uniEntry['Id2'] = re.split('\s+', currLine)[1]
				uniEntry['ReviewStatus'] = re.split('\s+', currLine)[2]
				uniEntry['length'] = " ".join(re.split('\s+', currLine)[3:]).split(".")[0]
			if (re.split('\s+', currLine)[0]=="AC" and ID1_empty):
				uniEntry['Id1'] = re.split('\s+', currLine)[1].split(';')[0]
				ID1_empty = 0
			if re.split('\s+', currLine)[0]=="OX":
				uniEntry['taxon'] = re.search('\d+', currLine).group(0)
			if re.split('\s+', currLine)[0]=="OS":
				species = " ".join([species, " ".join(re.split('\s+', currLine)[1:])])
				uniEntry['species'] = species[1:]
			if re.split('\s+', currLine)[0]=="OC":
				l = re.split('\s+', currLine)
				l.pop(0)
				lineage = " ".join([lineage, " ".join(l)])
				uniEntry['lineage'] = lineage[1:]
			if re.split('\s+', currLine)[0]=="PE":
				uniEntry['evidence'] = " ".join(re.split('\s+', currLine)[1:])
			if re.split('\s+', currLine)[0]=="DR" and re.search('Ensembl;', currLine) and re.search('G0', currLine):
				for id in re.split('\s+', currLine):
					if re.search('G0', id):
						uniEntry['EnsID'] = re.search('\w+\d+', id).group(0)
						break				
	f.close()

def printEntry(uniEntry):
	if re.search('Eutheria', uniEntry['lineage']):
		out = open('uniprot_sprot.parsed.eutheria2', 'a')
		out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(uniEntry['Id1'], uniEntry['Id2'], uniEntry['length'], uniEntry['taxon'], uniEntry['species'], uniEntry['lineage'], uniEntry['evidence'], uniEntry['ReviewStatus'], uniEntry['EnsID']))
		out.close()

if __name__ == "__main__":
	args = readArgs()
	run(**args)
