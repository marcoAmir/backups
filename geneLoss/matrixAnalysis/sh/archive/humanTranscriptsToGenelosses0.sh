#!/bin/bash -e

src=$GENELOSS/src
data=$GENELOSS/data

ensemblDir=/cluster/u/amirma/geneLoss/data/Ensembl
ucscDir=/cluster/u/amirma/geneLoss/data/UCSC
callTable=/cluster/u/amirma/geneLoss/matrixAnalysis/data/hg19.calltable
callTablePlacental=/cluster/u/amirma/geneLoss/matrixAnalysis/data/hg19.calltable.placentals
callTablePlacentalUnique=$data/hg19.calltable.unique.placentals
callVec=/cluster/u/amirma/geneLoss/matrixAnalysis/data/hg19.vec
callTableColumns=/cluster/u/amirma/geneLoss/matrixAnalysis/data/columns.hg19.calltable
speciesList=/cluster/u/amirma/geneLoss/matrixAnalysis/data/speciesList.txt

hfiltersApplied='h0'
########################################## FLAGS: ####################################################

segDupOverlap=0

transCompl=$1			# flag 1: filter out human transcript whose cdsStartStat or 
				# cdsEndStat are not complete ('cmpl')

transComplUtr=$2		# flag 2: filter out human transcript whose cdsStartStat or 
				# cdsEndStat are not complete ('cmpl'), and txStart=cdsStart or txEnd=cdsEnd

rmSegDup=$3			# flag 3: Whether to exclude from downstream analysis genes 
				# overlapping seg-duplication regions

mutualLoss=$4			# flag 4: Whether to to focus on gene losses in entire clades 
				# (i.e., losses shared >=2 neighboring taxa)

deduplicatedAnnotations=$5	# flag 5: Whether to remove all species-transcript annotations
				# where calls were made from multiple chains

permisiveLoF=$6			# flag 6: if selected, this flag will consider any of: framshift, stop-codon, 40% exon-deleted
				# as a transcript Lof (under the assumption that LoF could happen by small insertion for instance)
				# Current loosen default settings are: at least 2 framshift OR at least 2 stop-codons OR 40% deleted exons
				# To modify, refer to awk commands at lines 51 & 58

removeCurrentlyNonCoding=$7	# flag 7: There are about 511 protein-coding from ensembl-biomart74 that are no longer protein-coding in the latest ensembl version
				# if selected, this will exclude these 511 genes from the analysis 


########################################################################################################

if [ "$#" -ne 7 ]; then
  echo -e "\nUsage:\n  $0 h1 h2 h3 h4 h5 h6 h7\n\th1 - (choose 1 or 0) filter out transcripts with incomplete start/end\n\th1 - (choose 1 or 0) like h1 but in addition requires UTRs\n\th3 - (choose 1 or 0) removes transcripts overlapping segmental duplication regions\n\th4 - (choose 1 or 0) subset to gene losses shared by entire clades\n\th5 - (choose 1 or 0) whether to ignore all the species-transcript pairs with multiple annotation (i.e., calls from multiple chains)\n\th6 - (choose 1 or 0) if selected, a more permisive transcript loss-of-function applied, where any frameshift, premature stop codon or exon deletin will mark the transcript with Lof\n\th7 - Whether to exclude genes that are no longer 'protein-coding' under the most recent ensembl version"
  exit 1
fi

# (I) re-processs a ternary vector of transcripts using the flagged preferences:

if [ $deduplicatedAnnotations -ne 0 ]
then
	if [ $permisiveLoF -ne 0 ]; then
		awk -F'\t' '{if((NR>1) && ($4>3 || ($10/$14>0.33) || $11>3)) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t1"} else {print $0}}' $callTablePlacentalUnique > tmpLoose
		callTablePlacentalUnique=tmpLoose
	fi
	echo $callTablePlacentalUnique
	$src/sh/callsToVec.sh $callTablePlacentalUnique hg19.vec
else
	if [ $permisiveLoF -ne 0 ]; then
		awk -F'\t' '{if((NR>1) && ($4>3 || ($10/$14>0.33) || $11>3)) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t1"} else {print $0}}' $callTable > tmpLoose
		callTable=tmpLoose
	fi
	echo $callTable
	$src/sh/callsToVec.sh $callTable hg19.vec
fi

join -t$'\t' -1 1 -2 1 -o '2.1 2.2 2.3 2.4' <(awk -F'\t' '{if($6=="protein_coding") print $2}' $ensemblDir/geneTranscript.ensembl74.tsv | sort -u) <(sort hg19.vec) > tmp; mv tmp hg19.vec	# by default, I'll keep only transcripts belonging to protein-coding genes
																								# this removes 22 transcripts belonging to 'polymorphic_pseudogene'

if [ $transCompl -ne 0 ]
then
	#join -t$'\t' -1 1 -2 1 -o '2.1 2.2 2.3 2.4' <(hgsql hg19 -Ne "SELECT name FROM hg19.ensGene WHERE cdsStartStat='cmpl' AND cdsEndStat='cmpl'" | sort -u) <(sort hg19.vec) > tmp; mv tmp hg19.vec
	join -t$'\t' -1 1 -2 1 -o '2.1 2.2 2.3 2.4' <(sort -u $data/Ensembl/hg19.CompleteTranscripts) <(sort hg19.vec) > tmp; mv tmp hg19.vec
	hfiltersApplied=$hfiltersApplied,'h1'
fi

if [ $transComplUtr -ne 0 ]
then
	#join -t$'\t' -1 1 -2 1 -o '2.1 2.2 2.3 2.4' <(hgsql hg19 -Ne "SELECT name FROM hg19.ensGene WHERE cdsStartStat='cmpl' AND cdsEndStat='cmpl' AND cdsStart!=txStart AND cdsEnd!=txEnd" | sort -u) <(sort hg19.vec) > tmp; mv tmp hg19.vec
	join -t$'\t' -1 1 -2 1 -o '2.1 2.2 2.3 2.4' <(sort -u $data/Ensembl/hg19.CompleteTranscriptsUTRs) <(sort hg19.vec) > tmp; mv tmp hg19.vec
	hfiltersApplied=$hfiltersApplied,'h2'
fi

if [ $rmSegDup -ne 0 ]
then
	#hgsql hg19 -Ne "SELECT chrom, cdsStart, cdsEnd, name FROM hg19.ensGene" > a.bed
	#hgsql hg19 -Ne "SELECT chrom, chromStart, chromEnd, name FROM hg19.genomicSuperDups" > b.bed
	cp $data/Ensembl/hg19.transcripts.bed a.bed
	cp $data/UCSC/hg19.segmentalDuplications.bed b.bed
	# join -t$'\t' -1 1 -2 1 -o '2.1 2.2 2.3 2.4' -v 2 <(overlapSelect a.bed b.bed -statsOutput stdout | awk -F'\t' '{if($4>'$segDupOverlap') print $0}' | cut -f2 | sort -u) <(sort hg19.vec) > tmp; mv tmp hg19.vec
	join -t$'\t' -1 1 -2 1 -o '2.2' <(join -t$'\t' -1 1 -2 2 -o '2.1' <(overlapSelect a.bed b.bed -statsOutput stdout | awk -F'\t' '{if($4>'$segDupOverlap') print $0}' | cut -f2 | sort -u) <(sort -t$'\t' -k2,2 $data/Ensembl/geneTranscript.ensembl74.tsv) | sort -u) <(sort -t$'\t' $data/Ensembl/geneTranscript.ensembl74.tsv) | sort -u > transcriptsInSegDup.txt
	join -t$'\t' -1 1 -2 1 -v 2 -o '2.1 2.2 2.3 2.4' <(sort transcriptsInSegDup.txt) <(sort hg19.vec) > tmp; mv tmp hg19.vec
	rm -rf a.bed b.bed transcriptsInSegDup.txt
	hfiltersApplied=$hfiltersApplied,'h3'
fi


# (II) - convert transcripts to gene losses in species
cat <(cut -f1,2 hg19.vec | sed -e 's/\t/\,/g' | awk -F',' '{for(i=2; i<=NF; i++) {print $1"\t"$i"\t0"}}') <(cut -f1,3 hg19.vec | sed -e 's/\t/\,/g' | awk -F',' '{for(i=2; i<=NF; i++) {print $1"\t"$i"\t?"}}') <(cut -f1,4 hg19.vec | sed -e 's/\t/\,/g' | awk -F',' '{for(i=2; i<=NF; i++) {print $1"\t"$i"\t1"}}') | awk -F'\t' '{if($2!="") print $0}' | sort > geneTranscriptqSpeciesCall.tsv

join -t$'\t' -1 3 -2 1 -o '1.1 1.2 2.1 2.2 2.3' <(join -t$'\t' -1 2 -2 1 -o '1.1 2.1 2.2' <(sort -t$'\t' -k2 ${ensemblDir}/humanGeneSymbol.humanEnsembl.biomart74.NoSyn.map) <(sort ${ensemblDir}/geneTranscript.ensembl74.tsv) | sort -t$'\t' -k3) <(sort geneTranscriptqSpeciesCall.tsv) | awk -F'\t' '{if($3!=$4) print $0}' | sort -u | sort -k2,2 -k4,4 > tmp; mv tmp geneTranscriptqSpeciesCall.tsv

#join -t$'\t' -1 3 -2 1 -o '1.1 1.2 1.3 1.4 1.5 2.2 2.3' <(sort -t$'\t' -k3 geneTranscriptqSpeciesCall.tsv) <(hgsql hg19 -Ne "SELECT name, cdsStartStat, cdsEndStat FROM hg19.ensGene" | sort) | sort -k2,2 -k4,4 > tmp; mv tmp geneTranscriptqSpeciesCall.tsv 
join -t$'\t' -1 3 -2 1 -o '1.1 1.2 1.3 1.4 1.5 2.2 2.3' <(sort -t$'\t' -k3 geneTranscriptqSpeciesCall.tsv) <(sort $data/Ensembl/hg19.transcripts.cdsStartEndStat) | sort -k2,2 -k4,4 > tmp; mv tmp geneTranscriptqSpeciesCall.tsv

sed -e 's/?/0/g' geneTranscriptqSpeciesCall.tsv | awk -F'\t' '{if(a[1]==$2"_"$4) {printf "\t"$5} else {a[1]=$2"_"$4; printf "\n"a[1]"\t"$5}}' | tail -n+2 | awk -F'\t' '{a[1]=0; for(i=2; i<=NF; i++) {a[1]=a[1]+$i}; if(a[1]/(NF-1)==1) {print $1"\t"$2}; if(a[1]/(NF-1)<1) {print $1"\t0"}}' | sed -e "s/_/\t/g" | sort -u > geneSpecies.call	# Convert the latter table to a simplified gene(ensembl ID)-query species gene loss calls. All unknowns will be treated as intacts 0->?: 
			# (0=gene is active; 1=gene inactivated). For a gene-species to be considered inactivated (1), all its transcripts should be inactivated. 

if [[ $removeCurrentlyNonCoding -ne 0 ]]
then
	join -t$'\t' -1 1 -2 1 -o '1.1 1.2 1.3' -v 1 <(sort geneSpecies.call) <(sort $data/Ensembl/genes74noLongerProteinCoding.txt) | sort -u > tmp; mv tmp geneSpecies.call
fi

join -t$'\t' -1 2 -2 1 -o '1.1 1.2 1.3 2.1 2.2 2.3' <(sort -t$'\t' -k2 geneSpecies.call) <(sort -t$'\t' -k1 $speciesList) | egrep "placental|primate|rodent" | cut -f1-3 > geneSpecies.placentals.call

if [[ $mutualLoss -ne 0 ]]
then
	Rscript $src/R/carveCladeLosses.R
	cat <(sed -e 's/; /\t/g' geneCladeLosses.txt | awk -F'\t' '{for(i=2; i<=NF; i++) print $1"\t"$i"\t1"}' | sort -u) <(comm -3 -2 <(sort -u geneSpecies.placentals.call) <(sed -e 's/; /\t/g' geneCladeLosses.txt | awk -F'\t' '{for(i=2; i<=NF; i++) print $1"\t"$i"\t1"}' | sort -u) | awk -F'\t' '{print $1"\t"$2"\t0"}') | sort -u > tmp; mv tmp geneSpecies.placentals.call
	#join -t$'\t' -1 1 -2 1 -o '1.1 1.2 1.3' <(sort geneSpecies.placentals.call) <(sort geneCladeLosses.txt) | sort -u > tmp; mv tmp geneSpecies.placentals.call
	hfiltersApplied=$hfiltersApplied,'h4'
fi

if [[ $deduplicatedAnnotations -ne 0 ]]
then
	hfiltersApplied=$hfiltersApplied,'h5'
fi

if [[ $permisiveLoF -ne 0 ]]
then
	hfiltersApplied=$hfiltersApplied,'h6'
fi

if [[ $removeCurrentlyNonCoding -ne 0 ]]
then
	hfiltersApplied=$hfiltersApplied,'h7'
fi

# (III) - Run some counts, and update the filter log, and summary counts:

totalTranscripts=`cut -f3 geneTranscriptqSpeciesCall.tsv | sort -u | wc -l`
totalGenes=`cut -f2 geneTranscriptqSpeciesCall.tsv | sort -u | wc -l`
totalGeneLosses=`awk -F'\t' '{if($3==1) print $0}' geneSpecies.call | sort -u | wc -l`
placentalGeneLosses=`awk -F'\t' '{if($3==1) print $0}' geneSpecies.placentals.call | sort -u | wc -l`
nGenes=`awk -F'\t' '{if($3==1) print $1}' geneSpecies.placentals.call | sort -u | wc -l`
nSpecies=`awk -F'\t' '{if($3==1) print $2}' geneSpecies.placentals.call | sort -u | wc -l`

if [ ! -e transToGene.filtersCounts.log ]
then 
	echo -e "filters\ttotalTranscripts\ttotalGenes\tGeneLosses\tPlacentalLosses\tnGenes\tnSpecies" > transToGene.filtersCounts.log
fi

echo -e "$hfiltersApplied\t$totalTranscripts\t$totalGenes\t$totalGeneLosses\t$placentalGeneLosses\t$nGenes\t$nSpecies" >> transToGene.filtersCounts.log
sed -e 's/h0,//g' transToGene.filtersCounts.log > tmp; mv tmp transToGene.filtersCounts.log

rm -rf tmpLoose
