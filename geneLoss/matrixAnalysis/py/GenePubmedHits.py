#!/usr/bin/env python
#
# Count pubmed hits to Gene-symbols that are mapped to Ensembl Identifiers
# And output also up to 20 most recent pubmed ID hits
import argparse
import os
import re
import socket
import sys

from Bio import Entrez

def readArgs():
	parser = argparse.ArgumentParser('Add Pubmed query counts '
		'for a given list of gene-symbols')
	parser.add_argument('geneMapFile', metavar='humanGeneSymbol.humanEnsembl.biomart74.NoSyn.map',
		help='A table mapping gene symbols to Ensembl identifiers')
	parser.add_argument('-o', '--output', dest='output_file',
		help='Write output to file (default is to stdout).')
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	return args

def run(geneMapFile, output_file):
	g = open(geneMapFile)
	out = open(output_file, 'w')
	warnings = open('warnings.txt', 'w')
	Entrez.email = "amarcovitz@gmail.com"
	for line in g:
		print line
		gene = line.split()[0]
		handle = Entrez.esearch(db="pubmed", term=gene)
		record = Entrez.read(handle)
		if int(record["Count"]) < 50000:
			out.write("%s\t%d\t%s\n" %(gene, int(record["Count"]), '|'.join(record["IdList"])))
		else:
			print 'Warning! likely an ambigious gene name: '+gene
			out.write("%s\t%d\t%s\n" %(gene, int(record["Count"]), 'Warning! likely an ambigious gene name'))
	g.close()
	out.close()

if __name__ == "__main__":
	args = readArgs()
	run(**args)
	
	#handle = Entrez.esearch(db="pubmed", term="PITX1")
	#record = Entrez.read(handle)
	#print record
	#print record["Count"]
	#print record["IdList"]
	#print len(record["IdList"])
	#searchResults = pubmed.searchGeneAndTrait("PITX1")
