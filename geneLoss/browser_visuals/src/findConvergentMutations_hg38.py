#!/usr/bin/python2.7

# Author: Mike Gloudemans
# Date created: 10/9/2015

import pprint
import operator
import sys
import os
import optparse
import subprocess

import time

import copy

import BranchLengthScoring as BLS

import logging

import datetime

cluster_mode = True
convergentSoft = True	# Added - August 22

pp = pprint.PrettyPrinter(indent = 4)

paml_control_template = '''
    seqfile = convergence/control/{0}_{1}.aa * sequence data file name
    outfile = convergence/output/{0}_{1}.mp        * main result file
   treefile = convergence/control/{0}_{1}.trees  * tree structure file name

    seqtype = 2  * 0:nucleotides; 2:amino acids, 3:binary
      ncatG = 8  * # of categories in the dG model of rates
      nhomo = 0  * nonhomogeneous in calcualting P for branch
'''

### Parameters and constants
transcript_file = "/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/filteredTranscripts.hg38.txt"					# Added/edited by Amir
#transcript_file = "filteredTranscripts.txt" 
transcript_gene_mapping_file = "/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/nonDuplicatedCanonical.hg38.bed"			# Added/edited by Amir 
exon_file = "/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data/Ensembl/GRCh38/hg38.exonAlignments.bed.ref"				# Added/edited by Amir 
gene_symbol_file = "/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data/Ensembl/GRCh38/humanGeneSymbol.humanEnsembl.biomart86.NoSyn.map"	# Added/edited by Amir 

# The following species list contains 61 placental mammals with
# alignments to hg38. The species
# "monDom5","sarHar1","macEug2","ornAna1" were removed because
# they fall within the marsupials or beyond.

species_list = sorted(["hg38","ailMel1","bosTau8","calJac3","camFer1","canFam3","capHir1","cavPor3","cerSim1",
		"chiLan1","chlSab2","chrAsi1","conCri1","criGri1","dasNov3","echTel2","eleEdw1","eptFus1","equCab2",
		"eriEur2","felCat8","hetGla2","jacJac1","lepWed1","loxAfr3","macFas5","mesAur1","micOch1",
		"mm10","musFur1","myoDav1","myoLuc2","nomLeu3","ochPri3","octDeg1","odoRosDiv1","orcOrc1",
		"oryAfe1","oryCun2","otoGar3","oviAri3","panHod1","panTro5","papAnu2","ponAbe2","pteAle1","pteVam1","rheMac8",
		"rn6","saiBol1","sorAra2","speTri2","susScr3","triMan1","tupChi1","turTru2","vicPac2"])

# this is the list we use in hg19
#species_list_old = sorted(["hg19", "tupChi1","speTri2","jacJac1","micOch1","criGri1","mesAur1","mm10","rn5",
#		"hetGla2","cavPor3","chiLan1","octDeg1","oryCun2","ochPri3","susScr3","vicPac2",
#		"camFer1","turTru2","orcOrc1","panHod1","bosTau7", "capHir1",
#		"equCab2","cerSim1","felCat5","canFam3","musFur1","ailMel1","odoRosDiv1","lepWed1",
#		"pteAle1","pteVam1","myoDav1","myoLuc2","eptFus1","eriEur2","sorAra2","conCri1","loxAfr3",
#		"eleEdw1","triMan1","chrAsi1","echTel2","oryAfe1","dasNov3", "panTro4", "gorGor3", "ponAbe2", 
#		"nomLeu3", "rheMac3", "macFas5", "papHam1", "chlSab1", "calJac3",
#		"saiBol1", "otoGar3", "oviAri1"])

legal_amino_acids = "ACDEFGHIKLMNPQRSTVWY"


# Temporary blacklisting
'''
blacklist = ["musFur1", "ailMel1"]
for bl in blacklist:
	species_list.remove(bl)
'''

conservation_window_padding = 5		# Check this many bases on either side of the site to compute a conservation score
# TODO: Implement this cutoff if desired. Not yet used.
minimum_species_cutoff = 40		# Amino acid must be present in AT LEAST this many species for us to consider it

def main():
	target_groups, outgroups, window_conservation, position_conservation, excludeList = parse()

	### Part 1: Load non-sequence data
	transcript_list = load_transcripts_from_file(transcript_file)
	transcript_to_gene = load_gene_transcript_mapping_from_file(transcript_gene_mapping_file)
	exon_locations = load_exons_from_file(exon_file)
	gene_symbols = load_gene_symbols_from_file(gene_symbol_file)

	### Part 2: Find shared coding mutations, output results to file

	target_species = []
	for tg in target_groups:
		for ts in tg:
			target_species.append(ts)
			assert ts in species_list
	species_string = "_".join(target_species)

	geneSpeciesExcludes = {}										# Exclude gene-species update - July 6th
	convergentFileOutput = "/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/convergent_mutations_hg38/convergentMutations_{0}_{1}_{2}.txt".format(species_string, window_conservation, position_conservation)	# Exclude gene-species update - July 6th
	backgroundFileOutput = "/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/convergent_mutations_hg38/background/background_{0}_{1}_{2}.txt".format(species_string, window_conservation, position_conservation)	# Exclude gene-species update - July 6th
	divergentFileOutput = "/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/divergent_mutations_hg38/divergentMutations_{0}_{1}_{2}.txt".format(species_string, window_conservation, position_conservation)	# Exclude gene-species update - July 6th
	if convergentSoft:
		convergentFileOutput = "/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/convergent_mutations_hg38/convergentMutations_{0}_{1}_{2}.txt.convergentSoft".format(species_string, window_conservation, position_conservation)	# Exclude gene-species update - July 6th
		backgroundFileOutput = "/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/convergent_mutations_hg38/background/background_{0}_{1}_{2}.txt.convergentSoft".format(species_string, window_conservation, position_conservation)	# Exclude gene-species update - July 6th
		divergentFileOutput = "/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/divergent_mutations_hg38/divergentMutations_{0}_{1}_{2}.txt.convergentSoft".format(species_string, window_conservation, position_conservation)	# Exclude gene-species update - July 6th
		
	if not excludeList=='NONE':										# Exclude gene-species update - July 6th
		exLst = open(excludeList)									# Exclude gene-species update - July 6th
		for line in exLst:										# Exclude gene-species update - July 6th	
			if line.strip().split()[0] in geneSpeciesExcludes.keys():				# Exclude gene-species update - July 6th	
				geneSpeciesExcludes[line.strip().split()[0]].append(line.strip().split()[1])	# Exclude gene-species update - July 6th	
			else:											# Exclude gene-species update - July 6th	
				geneSpeciesExcludes[line.strip().split()[0]] = []				# Exclude gene-species update - July 6th	
				geneSpeciesExcludes[line.strip().split()[0]].append(line.strip().split()[1])	# Exclude gene-species update - July 6th	
		exLst.close()	# Exclude gene-species update - July 6th
		convergentFileOutput = "/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/convergent_mutations_hg38/convergentMutations_{0}_{1}_{2}_{3}.txt".format(species_string, window_conservation, position_conservation, 'exclude' + excludeList.split("_")[1].split(".")[0])	# Exclude gene-species update - July 6th
		backgroundFileOutput = "/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/convergent_mutations_hg38/background/background_{0}_{1}_{2}_{3}.txt".format(species_string, window_conservation, position_conservation, 'exclude' + excludeList.split("_")[1].split(".")[0])	# Exclude gene-species update - July 6th
		divergentFileOutput = "/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/divergent_mutations_hg38/divergentMutations_{0}_{1}_{2}_{3}.txt".format(species_string, window_conservation, position_conservation, 'exclude' + excludeList.split("_")[1].split(".")[0])	# Exclude gene-species update - July 6th

	with open(convergentFileOutput, "w") as w:	# Added/edited by Amir
		with open(backgroundFileOutput, "w") as wTestCounts:	# Added/edited by Amir
			# Write file header
			# TODO: Write a header to correspond w/ actual output. It is no longer correct!
#			w.write("\t".join(["geneID", "transcriptID", "gene_symbol", "chromosome", "codon_start", \
#					"codon_end", "exon_number", "human_aligned_position", "human_protein_length", "A0", \
#					 "nA0", "number_aligned", "fraction_conforming", "target_species", \
#						"target_amino_acids", "other_nonconforming_species", "other_amino_acids", "number_nonconforming_species", 
#						"conservation_score", "conservation_window_padding"]) + "\n")
			
			wd = open(divergentFileOutput, "w")

			# Loop through every transcript in the list, finding shared mutations.
			iter = -1
			
			total_positions = 0				# Total positions at which humans have an amino acid.
			total_dropped_by_insufficient_species = 0	
			total_dropped_by_missing_targets = 0
			total_dropped_by_conservation = 0
			total_dropped_by_aa_mismatch = 0
			total_dropped_by_pamp_issue = 0
			total_dropped_by_no_convergence = 0
			total_converged_positions = 0
			total_positions_tested = 0			# Total number of positions at which we applied convergence tests.

			for transcript in transcript_list:

				# If our transcript isn't in the transcript_to_gene list, it's
				# because we had multiple transcripts map to the same gene, so
				# just continue to the next transcript.
				if transcript not in transcript_to_gene:
					continue

				# Also get rid of transcripts if we don't know the corresponding
				# gene symbol.
				if transcript_to_gene[transcript] not in gene_symbols:
					continue

				iter += 1
				bin = int(transcript[-2:])
				#with open("../../output/protein_alignments/{0}/{1}.txt".format(bin, transcript)) as alignment_file:
				with open("/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/protein_alignments_hg38/{0}/{1}.txt".format(bin, transcript)) as alignment_file:	# Added/edited by Amir
					sequence_list = {}

					# First, get human amino acid sequence
					data = alignment_file.readline().strip().split()
					assert transcript == data[0]
					sequence_list["hg38"] = data[1]

					# Then get protein alignments for other animals
					species_seen = set([]) # If a species has already been seen, throw away that species for this transcript
							       # because it's not a 1-to-1 alignment.
					for line in alignment_file:
						data = line.strip().split()
						# skip to the next if you see a manually excluded transcript-species-chain (commented line: # at the begining)
						if line[0]=="#":
							continue

						assert transcript == data[1]

						species = data[0]

						# Discard species that are in the exclude list of the transcript	# Exclude gene-species update - July 6th
						if not excludeList=='NONE':						# Exclude gene-species update - July 6th
							if transcript in geneSpeciesExcludes.keys():			# Exclude gene-species update - July 6th
								if species in geneSpeciesExcludes[transcript]:		# Exclude gene-species update - July 6th
									continue					# Exclude gene-species update - July 6th

						# Don't mess around with species not on our list
						if species not in species_list:
							continue

						# Discard species appearing multiple times
						if species in species_seen:
							try:
								del sequence_list[species]
							except:
								pass
							continue
						species_seen.add(species)
						sequence_list[species] = convert_alignment_to_sequence(data[3])
						assert len(sequence_list[species]) == len(sequence_list["hg38"])

					if transcript not in exon_locations:
						continue
					codon_locations = get_codon_locations(exon_locations[transcript])

					# Number of codons in list should always be equivalent to number
					# of positions in amino acid sequence. If not, then we need to
					# take a look at the data to see what's up.
					if len(codon_locations) != len(sequence_list["hg38"]):
						pp.pprint(exon_locations[transcript])
						print(len(codon_locations))
						print()
						pp.pprint(exon_locations[transcript])
						pp.pprint(sequence_list["hg38"])
						print("Bad exon listing? Skipping to next sequence.")
						print()
						continue
						

					number_tested = 0 # Number of positions that we've tested in this gene.

					# Now loop through the sequence and analyze one position at a time.
					for i in range(len(sequence_list["hg38"])):

						divergent = 0	# Divergent substitution update - July 6th

						# We're not interested in stop codons at the moment.
						if sequence_list["hg38"][i] == "*":
							continue

						total_positions += 1

						# Human sequences should only contain amino acids. If not, quit and alert
						# because I'm not expecting anything else.
						if sequence_list["hg38"][i] not in legal_amino_acids:
							print("Error in sequence", sequence_list["hg38"], "at position", i)
							assert False
						
						# Make map of amino acids found in every species.
						# Delete species with invalid amino acids.
						species_amino_acids = {}
						for species in species_list:
							if species in sequence_list:
								if sequence_list[species][i] in legal_amino_acids:
									species_amino_acids[species] = sequence_list[species][i]

						# If not enough species aligned here, then just move on.
						total_aligned = len(species_amino_acids.items())					
						if total_aligned < minimum_species_cutoff:
							total_dropped_by_insufficient_species += 1
							continue

						# Update (Aug19): if convergentSoft = True, accept position if at least one species from each target is represented
						# Old version (before Aug18): Only count this position as testable if it's aligned in all of our
						# target species and at least one member of each outgroup.
						if not present_in_sufficient_species(species_amino_acids, target_groups, outgroups, convergentSoft):
							total_dropped_by_missing_targets += 1
							continue

						AA_counts = count_amino_acids(species_amino_acids)
						AA_sort = sorted(AA_counts.items(), reverse = True, key = operator.itemgetter(1))
						A0, nA0 = AA_sort[0][0], AA_sort[0][1]
						

						# Note: This pre-check for conservation is not necessarily rigorous, but I'm
						# using it right now for the "early round" in hopes that it will speed the analysis
						# process.
						'''
						if nA0 * 1.0 / (len(species_amino_acids.items())-len(target_species)) < position_conservation - 0.2:
							print "got 1"
							total_dropped_by_conservation += 1
							continue
						'''							

						BBLS_conservation = get_BBLS_conservation(species_amino_acids, A0, target_groups) 
						#BBLS_conservation = round(BBLS_conservation, 1)	# Added October 18th, 2016
						if BBLS_conservation < position_conservation:
							total_dropped_by_conservation += 1
							continue
						
						if window_conservation > 0.001:
							if get_conservation_score(sequence_list.values(), i, conservation_window_padding) < window_conservation:
								total_dropped_by_conservation += 1
								continue

						number_tested += 1
						total_positions_tested += 1

						# Check to see whether our target species have the same amino acid at this position.
						# If so, then proceed with analysis. If not, check divergent substitution.
						target_AA, same_AA = groups_have_same_amino_acid(target_groups, species_amino_acids, convergentSoft)	# Update - Aug22
						if not same_AA:
							total_dropped_by_aa_mismatch += 1
							divergent, divergentAA = checkDivergent(target_groups, outgroups, species_amino_acids, A0, convergentSoft)	# Divergent substitution update - July 6th
							if not divergent:								# Divergent substitution update - July 6th
								continue								# Divergent substitution update - July 6th	

						# If outgroups have identical amino acids to our target species, then there's no chance
						# of convergence, so don't bother with the PAML analysis in these cases. This should
						# save us a bunch of time.
						if not outgroups_have_different_amino_acids(target_AA[0], outgroups, species_amino_acids):
							total_dropped_by_no_convergence += 1
							continue			

						print(total_positions)

						# This is where the good stuff goes. Infer ancestral sequences and determine whether convergent evolution
						# has occurred.		

						aligned_species = species_amino_acids.keys()
						# TODO: Maybe port this code to a separate function
						# Output sequences to a file for analysis with PAML pamp
						species_indices = {}
						with open("/cluster/u/amirma/rot/mike/bin/paml4.8/convergence/control/{0}_{1}.aa".format(species_string, position_conservation) ,"w") as wConv:		# Added/edited by Amir
							wConv.write("{0} 1\n\n".format(len(species_amino_acids.keys())))
							index = 1
							for species in species_amino_acids.keys():
								wConv.write("{0}   {1}\n".format(species, species_amino_acids[species]))
								species_indices[species] = index
								index += 1
						# END TODO

						
						# Use tree_doctor to create the corresponding Newick tree for analysis with PAML pamp
						trim_newick_tree(species_amino_acids.keys(),species_string, position_conservation)
						run_paml_pamp(species_string, position_conservation)

						parent, ancestral_seqs, confidence = parse_pamp_results(total_aligned, species_string, position_conservation)
						if parent == None:	
							number_tested -= 1
							total_positions_tested -= 1
							total_dropped_by_pamp_issue += 1
							continue

						# Determine which is outgroup branchpoint for each target group.
						key_ancestors = []
						for index in range(len(target_groups)):
							group = target_groups[index] + outgroups[index]
							ancestor_lines = []
							for species in group:
								
								# We have to do this check because not all outgroup species
								# need to be represented.
								if species not in species_indices:
									continue

								current = species_indices[species]
								line = [current]
								while current in parent:
									line.append(parent[current])
									current = parent[current]
								ancestor_lines.append(line)

							for j in range(len(ancestor_lines[0])):
								present_in_all = True
								for al in ancestor_lines:
									if ancestor_lines[0][j] not in al:
										present_in_all = False

								if present_in_all:
									key_ancestors.append(ancestor_lines[0][j+1])
									break

						if not shows_convergence(key_ancestors, target_groups, ancestral_seqs, species_amino_acids, convergentSoft):
							total_dropped_by_no_convergence += 1
							continue
						
						if divergent==0:	# Divergent substitution update - July 6th
							total_converged_positions += 1
						# If we make it this far, convergent evolution seems to have occurred.

						# Make a list of all species besides the target species that also have substitutions in their sequences
						nonconforming_species = []
						nonconforming_AAs = []
						conforming_species = []
						for species in species_list:
							if species in sequence_list and sequence_list[species][i] != A0 and \
							 species not in target_species and sequence_list[species][i] not in "-?":
								nonconforming_species.append(species)
								nonconforming_AAs.append(sequence_list[species][i])
							else:
								conforming_species.append(species)
								

						# Make sure we have at least something to output in every case, so we don't break the file format.
						nonconforming_count = len(nonconforming_species)
						if nonconforming_count == 0:
							nonconforming_species = ["NA"]
							nonconforming_AAs = ["NA"]

						gene_symbol = gene_symbols[transcript_to_gene[transcript]]

						# TODO: Helper function?
						# Write data to file
						# TODO: Add header to clarify what the different columns mean.
						if not divergent:	# Divergent substitution update - July 6th
							w.write("\t".join([transcript_to_gene[transcript],
								  transcript,
								  gene_symbol,
								  codon_locations[i][0],
								  str(codon_locations[i][1]),
								  str(codon_locations[i][2]),
								  str(codon_locations[i][3]),
								  str(i),
								  str(len(sequence_list["hg38"])),
								  str(A0),
								  str(nA0),
								  str(total_aligned),
								  str(nA0 * 1.0 / total_aligned),
								  "{0}".format("|".join(target_species)),
								  target_AA[0],		# Older version:        "{0}".format(sequence_list[target_groups[0][0]][i]),
								  "|".join(nonconforming_species),
								  "|".join(nonconforming_AAs),
								  "|".join(conforming_species),
								  str(nonconforming_count),
								  str(BBLS_conservation),
#								  str(get_conservation_score(sequence_list.values(), i, conservation_window_padding)),
								  str(conservation_window_padding),
								  str(get_conservation_score(sequence_list.values(), i, 2)),
								  "2",
								  str(get_conservation_score(sequence_list.values(), i, 10)),
								  "10",
								  str(confidence)
								  ])
									+ "\n")
						else:
							wd.write("\t".join([transcript_to_gene[transcript],
								  transcript,
								  gene_symbol,
								  codon_locations[i][0],
								  str(codon_locations[i][1]),
								  str(codon_locations[i][2]),
								  str(codon_locations[i][3]),
								  str(i),
								  str(len(sequence_list["hg38"])),
								  str(A0),
								  str(nA0),
								  str(total_aligned),
								  str(nA0 * 1.0 / total_aligned),
								  "{0}".format("|".join(target_species)),
								  divergentAA,
								  "|".join(nonconforming_species),
								  "|".join(nonconforming_AAs),
								  "|".join(conforming_species),
								  str(nonconforming_count),
								  str(BBLS_conservation),
#								  str(get_conservation_score(sequence_list.values(), i, conservation_window_padding)),
								  str(conservation_window_padding),
								  str(get_conservation_score(sequence_list.values(), i, 2)),
								  "2",
								  str(get_conservation_score(sequence_list.values(), i, 10)),
								  "10",
								  str(confidence)
								  ])
									+ "\n")

					gene_symbol = gene_symbols[transcript_to_gene[transcript]]
					# Write data to file indicating how many amino acid positions
					# we actually tested in this gene, for use in later enrichment
					# analysis.
					wTestCounts.write("\t".join([transcript_to_gene[transcript],
						  transcript,
						  gene_symbol,
						  str(number_tested)]) + "\n")
	
	wd.close()	# Divergent substitution update - July 6th	
	os.chdir("/cluster/u/amirma/rot/mike/output_hg38")	# Added/edited by Amir
	with open("/cluster/u/amirma/git/forwardGenomics/geneLoss/convergentEvolution/convergent_mutations_hg38/summary/convergence_finder_summary_{0}_{1}_{2}.tmp".format(species_string,window_conservation,position_conservation), "w") as wSummary:	# Added/edited by Amir
		wSummary.write("Total positions: {0}\n".format(total_positions))
		wSummary.write("Total dropped by insufficient species count: {0}\n".format(total_dropped_by_insufficient_species))
		wSummary.write("Total dropped because missing in targets: {0}\n".format(total_dropped_by_missing_targets))
		wSummary.write("Total dropped by low conservation: {0}\n".format(total_dropped_by_conservation))
		wSummary.write("Total dropped by amino acid mismatch in targets: {0}\n".format(total_dropped_by_aa_mismatch))
		wSummary.write("Total dropped due to PAML pamp issue: {0}\n".format(total_dropped_by_pamp_issue))
		wSummary.write("Total dropped due to non-convergence: {0}\n".format(total_dropped_by_no_convergence))
		wSummary.write("Total positions converged: {0}\n".format(total_converged_positions))
		wSummary.write("Total positions tested: {0}\n".format(total_positions_tested))

# Function: convert_alignment_to_sequence
#
# Input: A protein sequence in aligned notation
#	e.g. "????A0b2C1D1D0D0D0---*"
#
# Output: The same protein sequence, with numbers removed
#       and all letters capitalized.
#	e.g. "????ABCDDDD---*"

def convert_alignment_to_sequence(alignment):
	sequence = ""
	for l in alignment:
		if l not in "0123":
			sequence += l

	return sequence.upper()

# Function: get_codon_locations
#
# Input: A list of tuples specifying exons in a transcript.
#	[ (chr, start, end, strand) , ... ]
#
# Exons must appear in the order that they appear in the protein;
# exons on the negative strand, therefore, appear in the reverse order
# of exons on the positive strand. 
#
# Output: A list of tuples specifying the location of each individual
# 	amino acid in the sequence, in order of their appearance in the
#	protein.
#	[ (chr, start, end, exon_number), ... ]
#

# NOTE that "codon_start" and "codon_end" refer to the left and right sides as shown
# in the genome browser, not the actual start and end of the 3 DNA nucleotides in the codon.
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
	
# Function: get_conservation_score
#
# Inputs: sequences = list of sequences, position = index of central position,
# 	padding = number of bases to examine on either side of the central position.
#
#
# Output: A score ranging from 0 to 1, showing the level of conservation.
#	For now, the conservation score is simply computed in the following manner:
#		- Compute consensus sequence. 
#		- Examine each position in each sequence that is within $padding pairs
#		  of the central position to see whether the position matches the 
#		  consensus sequence. Do not examine the central position itself.
#		- Add 1 to our score count if the position matches, and 0
#		  if it does not match, or if there is a ? or gap in the 
#		  consensus sequence.
#		- Divide our score by the total number of positions examined
#		  to compute our final score.
#
# TODO: Experiment with different ways of setting up this scoring function
# and see what gives the best results.
# One possibility would be not to penalize ? or -, but to just treat them
# as missing data and not factor them into the score.
#

# Get a conservation score based on the Bayesian
# Branch Length scoring method. This method is hopefully
# more agnostic to the topology of the tree.
def get_BBLS_conservation(species_amino_acids, A0, target_groups):

	# Remove target species from list so as not to bias results
	trimmed_AA_index = copy.deepcopy(species_amino_acids)
	for tg in target_groups:
		for ts in tg:
			if ts in trimmed_AA_index:
				del trimmed_AA_index[ts]

	species_to_keep = []
	for species in trimmed_AA_index.keys():
		species_to_keep.append(species)


	# Load and prune tree

	text = os.popen("tree_doctor /cluster/u/amirma/rot/mike/data/phylogeny/newick/mammals_hg38.nh -P {0}".format(",".join(species_to_keep))).read()	# Added/edited by Amir
	tree = BLS.parse_tree(text)

	'''
	with open("/cluster/u/amirma/rot/mike/bin/paml4.8/convergence.trees") as f:
		for i in range(2):
			f.readline()
		text = f.read().strip()
		tree = BLS.parse_tree(text)



	os.system(mammalian_tree_command)

	cmd = "tree_doctor pipe.tmp -n -P {0}".format(",".join(species_to_keep))

	p1 = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	text = p1.communicate()[0]

		
	tree = BLS.parse_tree(text)
	'''
	# Construct dictionary showing which species have the
	# canonical amino acid
	leaves = {}
	for key in trimmed_AA_index.keys():
		if trimmed_AA_index[key] == A0:
			leaves[key] = 1.0
		else:
			leaves[key] = 0.0

	# Compute BBLS and max possible BBLS
	bbls = BLS.BBLS(tree, leaves)
	max_bbls = BLS.getMaxBBLS(tree)

	return bbls / max_bbls


def get_conservation_score(sequences, position, padding):
	# Sanity check
	for seq in sequences:
		assert len(seq) == len(sequences[0])
	
	# Trim off the parts of the sequences that we actually care about
	cut_seqs = []
	for seq in sequences:
		# The extra math here ensures that we don't go beyond the edges of the sequence.
		cs = seq[max(0,position-padding):position] + seq[position+1:min(position+padding+1, len(sequences[0]))]
		cut_seqs.append(cs)

	consensus = compute_consensus_sequence(cut_seqs)

	# Compute score
	score = 0
	max_score = 0
	for i in range(len(cut_seqs[0])):
		for cs in cut_seqs:
			max_score += 1
			if cs[i] == consensus[i] and cs[i] not in "?-":
				score += 1
	return score*1.0/max_score

# Function: compute_consensus_sequence
#
# Input: a list of sequences of the same length
#
# Output: a single sequence showing the character that appears
#	most frequently at each position within the input
#	sequences.
#

def compute_consensus_sequence(sequences):
	# Sanity check
	for seq in sequences:
		assert len(seq) == len(sequences[0])
	
	# Build consensus sequence
	consensus = ""
	for i in range(len(sequences[0])):
		dCounts = {}
		for seq in sequences:
			dCounts[seq[i]] = dCounts.get(seq[i], 0) + 1
		consensus += sorted(dCounts.items(), reverse = True, key = operator.itemgetter(1))[0][0]

	return consensus

# Function: count_amino_acids
# 
# Count how many times each amino acid appears at this position in the sequence.
# Returns a dictionary of counts.
#
# Input: species_amino_acids = a dictionary showing which amino
#   acid is found in each aligning species
#
# Output: a dictionary of counts for each amino acid.
#
def count_amino_acids(species_amino_acids):
	AA_counts = {}
	for species in species_amino_acids.keys():
		AA = species_amino_acids[species]
		if AA not in AA_counts:
			AA_counts[AA] = 0
		AA_counts[AA] += 1
	return AA_counts	

# Function: load_exons_from_file
#
# Input: File containing locations of all exons
#   for each transcripts
#
# Output: Dictionary showing the genomic
#   locations of each exon within each transcript.
#
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

# Function: load_transcripts_from_file
#
# Input: File containing list of all transcripts
#   we want to use
#
# Output: List containing all transcripts.
#
def load_transcripts_from_file(transcript_file):
	transcript_list = []
	with open(transcript_file) as f:
		for line in f:
			if line[0]=="#":
				continue
			transcript_list.append(line.strip())
	return transcript_list

# Function: load_gene_transcript_mapping_from_file
#
# Input: File containing mapping of Ensembl transcripts to Ensembl
#   genes
#
# Output: Dict of mappings from Ensembl transcripts to Ensembl gene IDs
#
def load_gene_transcript_mapping_from_file(transcript_gene_mapping_file):
	transcript_to_gene = {}
	
	genes_mapped = set([])

	with open(transcript_gene_mapping_file) as f:
		for line in f:
			data = line.strip().split()
			# If two or more transcripts map to the same gene,
			# just toss them for now. This shouldn't be happening.
			if data[3] in genes_mapped:
				if data[4] in transcript_to_gene:
					del genes_mapped[transcript_to_gene]			
			else:
				transcript_to_gene[data[4]] = data[3]
			genes_mapped.add(data[3])
	return transcript_to_gene

# Function: load_gene_symbols_from_file
#
# Input: File containing mapping of genes to gene symbols
#
# Output: Dict mapping Ensembl gene names to gene symbols
#
def load_gene_symbols_from_file(gene_symbol_file):
	gene_symbols = {}
	with open(gene_symbol_file) as f:
		for line in f:
			data = line.strip().split()
			gene_symbols[data[1]] = data[0]
	return gene_symbols

# Function: groups_have_same_amino_acid
# NEW VERSION: 08/19/2016
# convergentSoft=True
# Can tolerate a missing sequence in some of the target species as long as:
# 	1) the 'convergent' amino acid is present at least once in all the targets
#	2) no other amino acids are present in the target species
#
# OLDER VERSION: prior to 08/18/2016
# (and if convergentSoft=False) Determine whether
#	1) amino acid is present in all species of target groups
#	2) amino acid is the same in all specis of target gropus
#
# Input: A list of target groups, each of which is a list
#   of target species of interest.
#
# Output: True if above criteria are satisified, otherwise False.
#
def groups_have_same_amino_acid(target_groups, species_amino_acids, convergentSoft):
	if convergentSoft:
		group1_AA=[]
		group2_AA=[]
		for species in target_groups[0]:
			if species in species_amino_acids:
				group1_AA.append(species_amino_acids[species])
		for species in target_groups[1]:
			if species in species_amino_acids:
				group2_AA.append(species_amino_acids[species])
		return list(set(group1_AA)), (list(set(group2_AA))==list(set(group1_AA)) and len(list(set(group1_AA)))==1 and len(list(set(group2_AA)))==1)
	else:
#	THE OLDER VERSION IS FROM HERE AND ON:
		for group in target_groups:
			for species in group:
				if species not in species_amino_acids or \
					species_amino_acids[species] != species_amino_acids[target_groups[0][0]]:
					return [''], False
		return species_amino_acids[target_groups[0][0]], True


#
# Function: outgroups_have_different_amino_acids
#
# Input: amino acid found in target species, list of outgroups,
#   dict of amino acids at this position in each species.
#
# Output: True if for every outgroup, there exists a species
#   that does not have the target amino acid. False otherwise.
#
def outgroups_have_different_amino_acids(target_amino_acid, outgroups, species_amino_acids):
	for og in outgroups:
		success = False
		for species in og:
			if species in species_amino_acids and species_amino_acids[species] != target_amino_acid:
				success = True
				break
		if not success:
			return False

	return True

# Function: trim_newick_tree
#
# Call tree doctor and prune the mammalian Newick tree.
#
# Input: list of species that should be preserved
#   in the outputted newick tree.
#
def trim_newick_tree(species_to_keep, species_string, position_conservation):
	os.system("printf '{0}\\t1\\n\\n' > /cluster/u/amirma/rot/mike/bin/paml4.8/convergence/control/{1}_{2}.trees".format(len(species_to_keep),species_string, position_conservation))	# Added/edited by Amir
	os.system("tree_doctor /cluster/u/amirma/rot/mike/data/phylogeny/newick/mammals_hg38.nh -P {0} >> /cluster/u/amirma/rot/mike/bin/paml4.8/convergence/control/{1}_{2}.trees".format(",".join(species_to_keep), species_string, position_conservation))	# Added/edited by Amir

# Function: run_paml_pamp
#
# Run PAML pamp
def run_paml_pamp(species_string, position_conservation):
	
	# Write control file
	with open("/cluster/u/amirma/rot/mike/bin/paml4.8/convergence/control/{0}_{1}.ctl".format(species_string, position_conservation), "w") as w:	# Added/edited by Amir
		w.write(paml_control_template.format(species_string, position_conservation))
	# Run PAML

	os.chdir("/cluster/u/amirma/rot/mike/bin/paml4.8")	# Added/edited by Amir
	run_command("/cluster/u/amirma/rot/mike/bin/paml4.8/pamp /cluster/u/amirma/rot/mike/bin/paml4.8/convergence/control/{0}_{1}.ctl".format(species_string, position_conservation))	# Added/edited by Amir
	os.chdir("/cluster/u/amirma/rot/mike/scripts/convergence")	# Added/edited by Amir


# Function: parse_pamp_results
#
# Read results of a PAML pamp run and extract
#   critical information.
#
# Output: parent = dict representation of phylogenetic
#   tree. ancestral_seqs = amino acids inferred at each
#   ancestral node.
#
def parse_pamp_results(total_aligned, species_string, position_conservation):
	parent = {}
	with open("/cluster/u/amirma/rot/mike/bin/paml4.8/convergence/output/{0}_{1}.mp".format(species_string, position_conservation)) as f:	# Added/edited by Amir
		# First figure out tree structure from file
		line = f.readline()
		while line != "":
			if not line.startswith("(1) Branch lengths and substitution pattern"):
				line = f.readline()
				continue
			break
		tree = f.readline().strip().split()
		for t in tree:
			data = [int(d) for d in t.split("..")]
			parent[data[1]] = data[0]	


		while line != "":
			if not line.startswith("(3) Parsimony reconstructions"):
				line = f.readline()
				continue
			break
		for i in range(4):
			f.readline()
		line = f.readline()
		seqs = line.split("|")[-1].strip()
		if seqs == "":
			seqs = line.split("|")[0].split(":")[1].strip()
		try:
			confidence = float(seqs.split("(")[1].split(")")[0])
			print(confidence)
			seqs = seqs.split("(")[0].strip()
		except:
			confidence = 1.0

		ancestral_seqs = {}
		issue_detected = False
		for j in range(len(seqs)):
			if seqs[j] not in legal_amino_acids:
				issue_detected = True
				print(seqs[j])
				break
			ancestral_seqs[total_aligned + j + 1] = seqs[j]

		if issue_detected:
			# Throw an error if we have an issue; in that case
			# we'll just move on the next sequence.
			return (None, None, None)

	return (parent, ancestral_seqs, confidence)

# Function: shows_convergence
#
# Input: list of indices of key ancestors, list of species groups in which we
#   are testing for convergence, list of ancestral sequences at internal
#   tree nodes, list of amino acids at each extant species.
#
# Output: True if the target groups have converged,
#   otherwise False.

def shows_convergence(key_ancestors, target_groups, ancestral_seqs, species_amino_acids, convergentSoft):
	# Update (Aug19):
	# We just want to make sure that the ancestral amino acid is not present in any of the target groups, but we tolerate missing sequence
	if convergentSoft:
		for i in range(len(key_ancestors)):
			group_AA=[]
			for species in target_groups[i]:
				if species in species_amino_acids:
					group_AA.append(species_amino_acids[species])
			if ancestral_seqs[key_ancestors[i]] in group_AA:
				return False
		return True
	# OLDER VESRION (Before Aug18):	
	else:
		for i in range(len(key_ancestors)):
			for tg in target_groups[i]:
				if tg in species_amino_acids:
					if ancestral_seqs[key_ancestors[i]] == species_amino_acids[target_groups[i][0]]:
						return False
				else:
					return False
		return True

# Function: present_in_sufficient_species
#
# Input: Amino acids at all aligning species, target species groups,
#   list of species in each outgroup.
#
# Output: (update Aug19)
#   if convergentSoft=True, retrun True if amino acid is aligned in at least one species of each target group 
# Older version output (Before Aug18), or if convergentSoft=False:
#   True if amino acid is aligned in all target species and
#   in at least one member of each outgroup, otherwise False.
#
def present_in_sufficient_species(species_amino_acids, target_groups, outgroups, convergentSoft):
	# New version: Make sure at least one member of each target_group is present.
	if convergentSoft:
		for tg in target_groups:
			present = False
			for ts in tg:
				if ts in species_amino_acids:
					present = True
					break
			if not present:
				return False
	else:
	# Make sure all target species are present.
		for tg in target_groups:
			for ts in tg:
				if ts not in species_amino_acids:
					return False

	# Make sure at least one member of each outgroup is present.
	for og in outgroups:
		present = False
		for os in og:
			if os in species_amino_acids:
				present = True
				break
		if not present:
			return False
		
	return True

# Function: deep_sort
#
# Input: A list of lists
#
# Output: Sorts the lower level of lists, and then sorts the
#   higher level. Useful for standardizing results.
#
#
def deep_sort(targets, outgroups):

	new_list = []

	for i in range(len(targets)):
		target_group = sorted(targets[i])
		outgroup = sorted(outgroups[i])
		new_list.append((target_group,outgroup))

	new_list = sorted(new_list)

	return [nl[0] for nl in new_list], [nl[1] for nl in new_list]

# Function: parse
#
# Parse command line arguments.
#
def parse():
	# Read command-line arguments
        parser = optparse.OptionParser(description="Detect convergent mutations in target species")
        parser.add_option('--tg', dest="target_groups", action="store")
        parser.add_option('--og', dest="outgroups", action="store")
        parser.add_option('--wc', dest="window_conservation", type=float, action="store", default=0.0)
        parser.add_option("--pc", dest="position_conservation", type=float, action="store", default=0.0)
	parser.add_option('--el', dest="excludeList", action="store")	# Exclude gene-species update - July 6th
        args = parser.parse_args()

        target_groups = []
        data = args[0].target_groups.strip().split(":")
        for d in data:
                target_groups.append(d.split(","))

        outgroups = [] 
        data = args[0].outgroups.strip().split(":")
        for d in data:
                outgroups.append(d.split(","))

        target_groups, outgroups = deep_sort(target_groups, outgroups)

        for tg in target_groups:
                for ts in tg:
                        assert ts in species_list
        for og in outgroups:
                for os in og:
                        assert os in species_list

	return (target_groups,outgroups,args[0].window_conservation, args[0].position_conservation, args[0].excludeList)	# Exclude gene-species update - July 6th

def run_command(command):
	subprocess.check_call(command, shell=True)

def checkDivergent(target_groups, outgroups, species_amino_acids, A0, convergentSoft):		# Divergent substitution update - July 6th
	# This subroutine makes sure that amino acids are different between the two target groups
	# and that the amino acid of each target is different than its outgroup
	#t1_agreement = True
	#t2_agreement = True
	#for s in target_groups[0]:
	#	t1_agreement = t1_agreement and (species_amino_acids[s]==species_amino_acids[target_groups[0][0]])
	#for s in target_groups[1]:
	#	t2_agreement = t2_agreement and (species_amino_acids[s]==species_amino_acids[target_groups[1][0]])
	
	#if not (t1_agreement and t2_agreement):
	#	return 0, 'X'
	aAcids1 = [];	aAcids2 = [];	aAcids1_og = [];	aAcids2_og = [];
	missingTargetAlignments = False	
	for s in target_groups[0]:
		if s in species_amino_acids:
			aAcids1.append(species_amino_acids[s])
		else:
			missingTargetAlignments = True
	for s in target_groups[1]:
		if s in species_amino_acids:
			aAcids2.append(species_amino_acids[s])
		else:
			missingTargetAlignments = True
	for s in outgroups[0]:
		if (s in species_amino_acids.keys()):
			aAcids1_og.append(species_amino_acids[s])
	for s in outgroups[1]:
		if (s in species_amino_acids.keys()):
			aAcids2_og.append(species_amino_acids[s])
	# Now, check that each outgroup is different than its target
	success_1 = (len(list(set(aAcids1) & set(aAcids1_og)))==0)
	success_2 = (len(list(set(aAcids2) & set(aAcids2_og)))==0)
	#if (success_1 and success_2 and len(list(set(aAcids1) & set(aAcids2)))==0):
	if convergentSoft:
		if (len(list(set(aAcids1) & set(aAcids2)))==0 and len(list(set(aAcids1) & set(A0)))==0 and len(list(set(aAcids2) & set(A0)))==0 and success_1 and success_2 and len(aAcids1)>0 and len(aAcids2)>0):
			return 1, '|'.join(aAcids1 + aAcids2)
		else:
			return 0, 'X'
	else:
		if (len(list(set(aAcids1) & set(aAcids2)))==0 and len(list(set(aAcids1) & set(A0)))==0 and len(list(set(aAcids2) & set(A0)))==0 and success_1 and success_2 and len(aAcids1)>0 and len(aAcids2)>0 and not missingTargetAlignments):
			return 1, '|'.join(aAcids1 + aAcids2)
		else:
			return 0, 'X'

if __name__ == "__main__":

	if cluster_mode:
		LOG_FILENAME = 'tmp/{}.txt'.format(datetime.datetime.now().strftime("%Y%m%d%H%M%S%f"))
		logging.basicConfig(filename=LOG_FILENAME,level=logging.DEBUG,)
		try:
			main()
		except:
			logging.exception('Got exception on main handler')
   			raise
	else:
		main()
