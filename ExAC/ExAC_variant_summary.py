#!/usr/bin/python

# Following documentation at:
#	http://pyvcf.readthedocs.io/en/latest/INTRO.html
# This script takes a subset of ExAC variants and outputs a summary for them

import vcf
import sys
import pandas as pd


populations = ['AFR', 'SAS', 'AMR', 'FIN',
	'EAS', 'NFE', 'OTH']

class Variant:
	def __init__(self, record):
		self.population = {}
		self.pos = str(record.CHROM) + ":" + str(record.POS)
		self.ref = record.REF
		self.alt = record.ALT
		self.info = record.INFO
		self.gene = {}
		for pop in populations:
			self.population[pop] = [int(self.info['AC_' + pop][0]), 
						int(self.info['AN_' + pop][0]),
		 				int(self.info['Hom_' + pop][0]), 
						float(self.info['AC_' + pop][0])/int(self.info['AN_' + pop][0])]
		self.population['Total'] = [int(self.info['AC_Adj'][0]), 
					int(self.info['AN_Adj'][0]),
		 			int(self.info['AC_Hom'][0]), 
					float(self.info['AC_Adj'][0])/int(self.info['AN_Adj'][0])]
		for var_type in self.info['CSQ']:
			var_type = var_type.split('|')
			self.gene[var_type[6]] = [var_type[3], var_type[4], 
			 var_type[7], var_type[1], var_type[2]]
	def gene_info(self, exclude_non_coding_transcripts=True):
		g = []
		for k,v in self.gene.items():
			if exclude_non_coding_transcripts:
				if v[2]=='protein_coding':
					g.append(k + "|" + v[0] + "|" + v[2] + "|" + v[3])
			else:
				g.append(k + "|" + v[0] + "|" + v[2] + "|" + v[3])
		return ";".join(g)
	def pop_info(self):
		populations_with_alt = []
		for pop in populations:
			if self.population[pop][0]>0:
				populations_with_alt.append(pop)
		return str(self.population['Total'][0]) + "\t" + str(self.population['Total'][1]) + "\t" + \
		 	str(self.population['Total'][2]) + "\t" + '{:0.3e}'.format(self.population['Total'][3]) + \
		 	"\t" + ";".join(populations_with_alt)
	def summary(self, population_summary=True, gene_info=True):
		print "\n\tVariant at " + self.pos 
		print "\t\tref allele: " + self.ref + "\n\t\talt allele: " + str(self.alt[0]) + "\n"
		if population_summary:
			print "\nVariant - Population Summary:\n"
			print pd.DataFrame(self.population, index=['Allele Count', 
			 'Allele Number', 'num. Homozygotes', 'Allele Freq.']).transpose()
		if gene_info:
			print "\nVariant - Genes Impacted:\n"
			print pd.DataFrame(self.gene, index=['Gene',
			 'EnsemblID', 'biotype', 'Variant Effect', 'Consequence']).transpose()


def iterate_vcf(in_vcf, chrom, start=None, end=None):
	vcf_reader = vcf.Reader(open(in_vcf, 'r')).fetch(chrom, start, end)
	for record in vcf_reader:
		curr_variant=Variant(record)
		pos = curr_variant.pos.split(":")
		print "chr" + pos[0] + "\t" + pos[1] + "\t" + curr_variant.ref + "\t" + \
		 str(curr_variant.alt) + "\t" + curr_variant.pop_info() + "\t" + curr_variant.gene_info()
	

if __name__ == "__main__":
	in_vcf, chrom, start, end = sys.argv[1:]
	iterate_vcf(in_vcf, int(chrom), int(start), int(end))
