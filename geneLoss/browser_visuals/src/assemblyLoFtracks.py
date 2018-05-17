#!/usr/bin/env python

# Amir Marcovitz
# Date created: 2/28/2017

# mammalian gene-losses: visualize mutations in the genome browser
# this script iterates over our list of transcripts and for a given target assembly generates bed files for all collected:
#	frameshift indels
#	stop-codons
#	non-synonymous mutation
#	deletions

import sys
import os
import argparse

import readDeletionsFromChains as RDEL

### Parameters and constants
transcript_file = "/cluster/u/amirma/geneLoss/hg38/browser_visuals/data/testedTranscripts.txt"
transcript_gene_mapping_file = "/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data/Ensembl/GRCh38/geneTranscript.ensembl86.tsv"
exon_file = "/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data/Ensembl/GRCh38/hg38.exonAlignments.bed.ref"
transcript_coding_exons = "/cluster/u/amirma/geneLoss/hg38/browser_visuals/data/testedTranscripts.ExonsCodingBases.bed"
gene_symbol_file = "/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data/Ensembl/GRCh38/humanGeneSymbol.humanEnsembl.biomart86.NoSyn.map"
chainMapDir = "/cluster/u/amirma/geneLoss/hg38/mappings/"

legal_amino_acids = "ACDEFGHIKLMNPQRSTVWY"
trancripts_indels = ''

def readArgs():
	print ""
	parser = argparse.ArgumentParser('Generate files for browser visualization of Loss of Function coding events')
	parser.add_argument('assembly', help='target assembly')
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	return args
	
def main(assembly):
	if not os.path.exists('./largeDel'):
    		os.makedirs('./largeDel')
	if not os.path.exists('./largeDel/' + assembly + '.gapsInChains.bed') or not os.path.exists('./largeDel/' + assembly + '.overlappingGaps'):
		subsetChains = True
		RDEL.readGapsInChainFile(assembly, subsetChains)
		RDEL.intersectDelWithCodingExons(assembly)
	print "\n\t...collecting mutations for target:  " + assembly + "\n"
	# Load non-sequence data
	transcript_list = load_transcripts_from_file(transcript_file)
	transcript_to_gene = load_gene_transcript_mapping_from_file(transcript_gene_mapping_file)
	exon_locations = load_exons_from_file(exon_file)
	gene_symbols = load_gene_symbols_from_file(gene_symbol_file)
	transcript_chainID = load_transcripts_chainMap(assembly)
	transcript_chainDel = get_transcript_chainDel(assembly, transcript_chainID)

	# Initiate output files:
	outDir = "/cluster/u/amirma/geneLoss/hg38/browser_visuals/browser_bedFiles/"
	nonSynFile = open(outDir + assembly + ".nonSyn.bed", 'w')
	stopCodonFile = open(outDir + assembly + ".stopCodon.bed", 'w')
	deletionsFile = open(outDir + assembly + ".deletions.bed", 'w')
	frameshiftsDelFile = open(outDir + assembly + ".del_frameshifts.bed", 'w')
	
	# Iterate over transcript list to extract all events
	print "\n\t...iterating over transcripts\n"
	tx_count = 0
	for transcript in transcript_list:
		tx_count += 1
		if transcript not in transcript_to_gene:
			print "  " + transcript + "not in transcript_to_gene\n"
			continue
		if transcript_to_gene[transcript] not in gene_symbols:
			print "  " + transcript + "not in gene_symbols\n"
			continue
		if transcript not in exon_locations:
			print "  " + transcript + "not in exon_locations\n"
			continue
		codon_locations = get_codon_locations(exon_locations[transcript])
		print str(tx_count) + "\t\t" + transcript
		bin = int(transcript[-2:])
		with open("/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/protein_alignments_hg38/{0}/{1}.txt".format(bin, transcript)) as alignment_file:
			sequence_list = {}
			# First, get human amino acid sequence
			data = alignment_file.readline().strip().split()
			assert transcript == data[0]
			sequence_list["hg38"] = data[1]
			for line in alignment_file:
				data = line.strip().split()
				# skip to the next if you see a manually excluded transcript-species-chain (commented line: # at the begining)
				if line[0]=="#":
					continue
				assert transcript == data[1]
				curr_assembly=data[0]
				if curr_assembly==assembly:
					trancripts_indels = data[2].split('|')
					sequence_list[assembly] = convert_alignment_to_sequence(data[3])
					assert len(sequence_list[assembly]) == len(sequence_list["hg38"])
		# Now loop through the sequence and analyze one position at a time.
		for i in range(len(sequence_list["hg38"])):
			# non-synonymous mutation:
			if sequence_list["hg38"][i]!=sequence_list[assembly][i] and sequence_list[assembly][i] in legal_amino_acids:
				mut = "".join([transcript, "_", sequence_list["hg38"][i], str(i), sequence_list[assembly][i]])
				nonSynFile.write("%s\t%d\t%d\t%s\n" %(codon_locations[i][0], codon_locations[i][1]-1, codon_locations[i][2], mut)) 
			if sequence_list[assembly][i]=="*" and sequence_list["hg38"][i] in legal_amino_acids:
				mut = "".join([transcript, "_", sequence_list["hg38"][i], str(i)])
				stopCodonFile.write("%s\t%d\t%d\t%s\n" %(codon_locations[i][0], codon_locations[i][1]-1, codon_locations[i][2], mut)) 
				# output codon to stop-codon bed file
				# print transcript + ":  stop-codon " + str(i) +  "  " + sequence_list["hg38"][i] + " " + sequence_list[assembly][i]
		if transcript in transcript_chainDel.keys():
			for curr_fragment in transcript_chainDel[transcript]:
				deletionsFile.write(curr_fragment + "\n")
		# I now want to figure out which deletion is a potential frameshift deletion
		if transcript in transcript_chainDel.keys():
			for i in range(len(trancripts_indels)):
				if "," in trancripts_indels[i]:
					indel_fragments = get_fs_dels(trancripts_indels[i], exon_locations[transcript][i], transcript_chainDel[transcript])
					for curr_fragment in indel_fragments:
						frameshiftsDelFile.write(curr_fragment + "\n")					
		# TO DO: show insertions
	nonSynFile.close()
	stopCodonFile.close()
	deletionsFile.close()
	frameshiftsDelFile.close()

def get_fs_dels(indels, exon_loc, chainDelFrags):
	# I currently ignore all question marks
	indel_fragments = []
	exon = list(exon_loc)
	exon[1] = exon[1]+6;   exon[2] = exon[2]-6	# the exons are flanked by 6bp non-exonic regions
	strand = exon[3]
	startPos = exon[1]
	# Iterate over chain gaps to see which is overlapping the exon and how much they delete
	for curr_fragment in chainDelFrags:
		x1 = int(curr_fragment.split("\t")[1])
		if x1>=exon[1] and x1<=exon[2]:
			l1 = min(int(curr_fragment.split("\t")[2]), exon[2])-x1
			if (l1 % 3 > 0):
				indel_fragments.append(curr_fragment)
	return indel_fragments



#	if strand == "-":
#		startPos = exon_loc[2]
#	for currBlock in indels.split(','):
#		if '-' in currBlock:
#			deletion = int(currBlock.split('-')[1])
#			if (deletion>0):
#				frameshift = deletion % 3
#		else:
#			if strand == "+":
#				startPos = startPos + int(currBlock)
#			if strand == "-":
#				startPos = startPos - int(currBlock)
				
def convert_alignment_to_sequence(alignment):
	sequence = ""
	for l in alignment:
		if l not in "0123":
			sequence += l
	return sequence.upper()
	
def get_codon_locations(exons):
	codon_locations = []
	# Each exon tuple contains 7 bases of 5' padding and 6 bases
	# of 3' padding on the Watson strand, as displayed in the genome browser.
	if exons[0][3] == "+":
		# Set current exon and position at beginning of gene
		done = False
		current_exon = 0
		start_pos = exons[current_exon][1] + 7
		stop_pos = start_pos + 2
		# Loop through exons, computing codon positions.
		while True:
			# If we're off the end of the exon, move to the next one.
			if start_pos > exons[current_exon][2] - 6:
				offset = start_pos - (exons[current_exon][2] - 6) - 1
				while start_pos > (exons[current_exon][2] - 6):
					current_exon += 1
					# If there is no next exon, then we're done.
					# Return the list of codon locations.
					if current_exon >= len(exons):
						return codon_locations
					start_pos = exons[current_exon][1] + 7 + offset
					stop_pos = start_pos + 2
			# Add codon position to list
			codon_locations.append((exons[0][0], start_pos, stop_pos, current_exon+1))
			# Advance to next codon position
			start_pos += 3
			stop_pos += 3
	else:
		# Set current exon and position at beginning of gene
		done = False
		current_exon = 0
		stop_pos = exons[current_exon][2] - 6
		start_pos = stop_pos - 2
		# Loop through exons, computing codon positions.
		while True:
			# If we're off the end of the exon, move to the next one.
			if stop_pos < exons[current_exon][1] + 7:
				offset = (exons[current_exon][1] + 7) - stop_pos - 1
				while stop_pos < exons[current_exon][1] + 7:
					current_exon += 1
					# If there is no next exon, then we're done.
					# Return the list of codon locations.
					if current_exon >= len(exons):
						return codon_locations
					stop_pos = exons[current_exon][2] - 6 - offset
					start_pos = stop_pos - 2
			# Add codon position to list
			codon_locations.append((exons[0][0], start_pos, stop_pos, current_exon+1))
			# Advance to next codon position
			start_pos -= 3
			stop_pos -= 3
	
# Some input-files reading functions:
	
def load_transcripts_from_file(transcript_file):
	transcript_list = []
	with open(transcript_file) as f:
		for line in f:
			if line[0]=="#":
				continue
			transcript_list.append(line.strip())
	return transcript_list

def load_gene_transcript_mapping_from_file(transcript_gene_mapping_file, protein_coding=True):
	transcript_to_gene = {}
	with open(transcript_gene_mapping_file) as f:
		for line in f:
			data = line.strip().split()
			if protein_coding:
				if data[-1]=="protein_coding":
					transcript_to_gene[data[1]] = data[0]
			else:
				transcript_to_gene[data[1]] = data[0]
	return transcript_to_gene
	
def load_exons_from_file(exon_file):
	exon_locations = {}
	# Then, load file w/ all exon locations
	with open(exon_file) as f:
		for line in f:
			data = line.strip().split()
			chrom = data[0]
			start = int(data[1])
			end = int(data[2])
			transcript = data[3].split(".")[0]
			strand = data[5]
			if transcript not in exon_locations:
				exon_locations[transcript] = []
			# For some reason, some of the genes in the input file
			# have the same exons listed multiple times.
			# Get rid of these duplicate exons.
			if (chrom, start, end, strand) not in exon_locations[transcript]:
				exon_locations[transcript].append((chrom, start, end, strand))
	return exon_locations
	
def load_gene_symbols_from_file(gene_symbol_file):
	gene_symbols = {}
	with open(gene_symbol_file) as f:
		for line in f:
			data = line.strip().split()
			gene_symbols[data[1]] = data[0]
	return gene_symbols
	
def load_transcripts_chainMap(assembly):
	transcript_chain = {}
	with open(chainMapDir + assembly + '.chain_ids') as f:
		for line in f:
			transcript_chain[line.strip().split(' ')[2]] = line.strip().split(' ')[1]
	return transcript_chain
	
def get_transcript_chainDel(assembly, transcript_chainID):
	print("\n\t...reading chain-gap overlapping transcripts\n")
	transcript_chainDel = {}
	for k in transcript_chainID.keys():
		transcript_chainDel[k] = []
	with open("./largeDel/" + assembly + '.overlappingGaps') as f:
		for line in f:
			curr_line = line.strip().split()
			if curr_line[0] in transcript_chainID.keys() and transcript_chainID[curr_line[0]]==curr_line[-1].split("_")[0]:
				transcript_chainDel[curr_line[0]].append("\t".join(curr_line[1:]))
	return transcript_chainDel

if __name__ == "__main__":
	args = readArgs()
	main(**args)
