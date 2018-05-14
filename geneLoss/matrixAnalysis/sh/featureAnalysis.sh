#!/bin/bash -e 

ensemblDir=/cluster/u/amirma/geneLoss/data/Ensembl
ucscDir=/cluster/u/amirma/geneLoss/data/UCSC
callTable=/cluster/u/amirma/geneLoss/matrixAnalysis/data/hg19.calltable
callTableColumns=/cluster/u/amirma/geneLoss/matrixAnalysis/data/columns.hg19.calltable
src=/cluster/u/amirma/geneLoss/matrixAnalysis/src

callsVsEnsembl=/cluster/u/amirma/geneLoss/matrixAnalysis/callsVsEnsembl.tsv

assemblyAcc=$1						# this is the assembly of the query species, with complete info, i.e., the number of the assembly
assembly=`echo $assemblyAcc | sed -e 's/[0-9]//g'`	# this will eliminate the number
species=$2
class=$3
mkdir $species						# at this point you'll have to remove previous directories
mkdir $species/TP $species/FP

hfiltersApplied='h0'	# which filters have been applied on the human side, to convert transcripts to gene calls
ofiltersApplied='o0'
####################################################### Filters: ################################################################################################################################################################
# (1) - partioning gene-loss calls into True and Flase positives
# On the human side we require for the gene counts that the transcripts are pristine protein-coding, intact transcripts (cdsStart/End are complete)
# On the ortholog species I only count protein-coding transcripts (by Ensembl biomart) and then can apply additional criteria  

orthoCompl=$4		# filter 1: a filter requiring that in the ortholog species the transcript is 
			# counted as false-positive only if both cdsStart/End are complete ('cmpl')

orthoComplUtr=$5	# filter 2: a filter requiring that in the ortholog species the transcript is 
			# counted as false-positive if the transcript has UTRs and the two of the above apply - (will override the preferences from above)

fracBasesDeleted=$6	# filter 3: a filter that count as gene losses - if the flag selected default is set to 25% deleted bases; otherwise, fraction is set to 0
frac=0


#################################################################################################################################################################################################################################

echo -e "\n  ...partioning Fp/Tp gene-loss calls for $species..."
# no flags
#join -t$'\t' -1 1 -2 1 -o '1.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10' <(awk -F'\t' '{if($3==1) print $0}' geneSpecies.placentals.call | grep -w $assembly | cut -f1 | sort -u) <(sort $ensemblDir/OrthoTables74/human${species}.EnsemblOrtho.biomart74.map) | awk -F'\t' '{IGNORECASE=1; if($4~/protein_coding/ && $9~/protein_coding/ && $3==$8) print $0}' | cut -f1,2 | sort -u > $species/FP/fpGeneTranscripts.tsv
join -t$'\t' -1 1 -2 1 -o '1.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10' <(awk -F'\t' '{if($3==1) print $0}' geneSpecies.placentals.call | grep -w $assembly | cut -f1 | sort -u) <(sort $ensemblDir/OrthoTables74/human${species}.EnsemblOrtho.biomart74.map) | awk -F'\t' '{IGNORECASE=1; if($4~/protein_coding/ && $9~/protein_coding/) print $0}' | cut -f1,2 | sort -u > $species/FP/fpGeneTranscripts.tsv

if [ $orthoCompl -ne 0 ]
then
	join -t$'\t' -1 7 -2 1 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10' <(join -t$'\t' -1 1 -2 1 -o '1.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10' <(awk -F'\t' '{if($3==1) print $0}' geneSpecies.placentals.call | grep -w $assembly | cut -f1 | sort -u) <(sort $ensemblDir/OrthoTables74/human${species}.EnsemblOrtho.biomart74.map) | awk -F'\t' '{IGNORECASE=1; if($4~/protein_coding/ && $9~/protein_coding/ && $3==$8) print $0}' | sort -t$'\t' -k7,7) <(hgsql ${assemblyAcc} -Ne "SELECT name, txStart, txEnd, cdsStart, cdsEnd, cdsStartStat, cdsEndStat FROM ${assemblyAcc}.ensGene WHERE cdsStartStat='cmpl' AND cdsEndStat='cmpl'" | cut -f1 | sort -u) | cut -f1,2 | sort -u > $species/FP/fpGeneTranscripts.tsv
	ofiltersApplied=$ofiltersApplied,'o1'
fi

if [ $orthoComplUtr -ne 0 ]
then
	join -t$'\t' -1 7 -2 1 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10' <(join -t$'\t' -1 1 -2 1 -o '1.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10' <(awk -F'\t' '{if($3==1) print $0}' geneSpecies.placentals.call | grep -w $assembly | cut -f1 | sort -u) <(sort $ensemblDir/OrthoTables74/human${species}.EnsemblOrtho.biomart74.map) | awk -F'\t' '{IGNORECASE=1; if($4~/protein_coding/ && $9~/protein_coding/ && $3==$8) print $0}' | sort -t$'\t' -k7,7) <(hgsql ${assemblyAcc} -Ne "SELECT name, txStart, txEnd, cdsStart, cdsEnd, cdsStartStat, cdsEndStat FROM ${assemblyAcc}.ensGene WHERE cdsStartStat='cmpl' AND cdsEndStat='cmpl' AND cdsStart!=txStart AND cdsEnd!=txEnd" | cut -f1 | sort -u) | cut -f1,2 | sort -u > $species/FP/fpGeneTranscripts.tsv
	ofiltersApplied=$ofiltersApplied,'o2'
fi

join -t$'\t' -1 1 -2 1 -o '2.1 2.2' <(comm -3 -2 <(awk -F'\t' '{if($3==1) print $0}' geneSpecies.placentals.call | grep -w $assembly | cut -f1 | sort -u) <(cut -f1 $species/FP/fpGeneTranscripts.tsv | sort -u) | sort -u) <(sort $ensemblDir/geneTranscript.ensembl74.tsv) | sort -u > $species/TP/tpGeneTranscripts.tsv

join -t$'\t' -1 2 -2 2 -o '1.1 1.2 2.14 2.15' <(sort -t$'\t' -k2 $species/TP/tpGeneTranscripts.tsv) <(sort -t$'\t' -k2 $ucscDir/hg19.ensGene.tsv) > tmp; mv tmp $species/TP/tpGeneTranscripts.tsv	# add the ensembl-transcript cdsStart, cdsEnd status (from santa cruz)
join -t$'\t' -1 2 -2 2 -o '1.1 1.2 2.14 2.15' <(sort -t$'\t' -k2 $species/FP/fpGeneTranscripts.tsv) <(sort -t$'\t' -k2 $ucscDir/hg19.ensGene.tsv) > tmp; mv tmp $species/FP/fpGeneTranscripts.tsv	# add the ensembl-transcript cdsStart, cdsEnd status (from santa cruz)
TPincompl=`tail -n+2 $species/TP/tpGeneTranscripts.tsv | awk -F'\t' '{if($3!="cmpl" || $4!="cmpl") print $2}' | sort -u | wc -l`
FPincompl=`tail -n+2 $species/FP/fpGeneTranscripts.tsv | awk -F'\t' '{if($3!="cmpl" || $4!="cmpl") print $2}' | sort -u | wc -l`
echo -e "     TP: $TPincompl incomplete transcripts"
echo -e "     FP: $FPincompl incomplete transcripts"

# Expand with the features of the transcript in each annotation
echo -e "\n  ...reading ortholgous TP-genes transcript calls for $species from hg19.calltable..."
awk -F'\t' '{for(i=1; i<=NF; i++) printf $i"\t"}' $callTableColumns  | awk -F'\t' '{print "ensgID\tenstID\t"$1"\t"$3"\t"$4"\t"$11"\t"$5"\t"$6"\t"$15"\t"$14"\tfraction.deleted.exons\tfraction.deleted.bases\tCDSStart\tCDSEnd"}' > tmp
join -t$'\t' -1 2 -2 2 -o '1.1 1.2 2.1 2.3 2.4 2.11 2.5 2.6 2.15 2.14 1.3 1.4' <(sort -k2 $species/TP/tpGeneTranscripts.tsv) <(grep $assemblyAcc $callTable | sort -k2,2) | sort -u | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$7/$9"\t"$8/$10"\t"$11"\t"$12}' >> tmp; mv tmp $species/TP/tpGeneTranscripts.tsv
echo -e "\n  ...reading ortholgous FP-genes transcript calls for $species from hg19.calltable..."
awk -F'\t' '{for(i=1; i<=NF; i++) printf $i"\t"}' $callTableColumns  | awk -F'\t' '{print "ensgID\tenstID\t"$1"\t"$3"\t"$4"\t"$11"\t"$5"\t"$6"\t"$15"\t"$14"\tfraction.deleted.exons\tfraction.deleted.bases\tCDSStart\tCDSEnd"}' > tmp
join -t$'\t' -1 2 -2 2 -o '1.1 1.2 2.1 2.3 2.4 2.11 2.5 2.6 2.15 2.14 1.3 1.4' <(sort -k2 $species/FP/fpGeneTranscripts.tsv) <(grep $assemblyAcc $callTable | sort -k2,2) | sort -u | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$7/$9"\t"$8/$10"\t"$11"\t"$12}' >> tmp; mv tmp $species/FP/fpGeneTranscripts.tsv

if [ $fracBasesDeleted -ne 0 ]
then
	frac=0.25
	awk -F'\t' '{if($12>'$frac') print $0}' $species/TP/tpGeneTranscripts.tsv > tmp; mv tmp $species/TP/tpGeneTranscripts.tsv
	awk -F'\t' '{if($12>'$frac') print $0}' $species/FP/fpGeneTranscripts.tsv > tmp; mv tmp $species/FP/fpGeneTranscripts.tsv
fi

# Add the count to the table:
echo -e "\n  ...updating precision table..."
if [ ! -e $callsVsEnsembl ]
then 
	echo -e "species\tlossCalls\tensemblProtCodOrtho\tgroup\thumanFilters\torthoFilters" > $callsVsEnsembl
fi
TPs=`cut -f1 $species/TP/tpGeneTranscripts.tsv | tail -n+2 | sort -u | wc -l`
FPs=`cut -f1 $species/FP/fpGeneTranscripts.tsv | tail -n+2 | sort -u | wc -l`
echo -e "$species\t$((TPs+FPs))\t$FPs\t$class\t$hfiltersApplied\t$ofiltersApplied\t$frac" | sed -e 's/\t0,/\t/g' >> $callsVsEnsembl
cat <(head -1 $callsVsEnsembl) <(tail -n+2 $callsVsEnsembl | sort -u | sort -k5 -k4) > tmp; mv tmp $callsVsEnsembl
sed -e 's/\th0,/\t/g;s/\to0,/\t/g' $callsVsEnsembl | sed -e "s/^Pig/#Pig/g" > tmp; mv tmp $callsVsEnsembl

# Now extract from the TP and FP, what genes overlap (any) segmental duplication event:
echo -e "\n  ...overlapping with seg-dup regions..."
join -t$'\t' -1 1 -2 1 -o '2.3 2.4 2.5 1.1' <(cut -f1 $species/FP/fpGeneTranscripts.tsv | sort -u) <(sort $ensemblDir/geneTranscript.ensembl74.tsv) | sort -u | overlapSelect stdin $ucscDir/hg19.segmentalDuplications.bed -selectFmt=bed -statsOutput stdout | tail -n+2 | cut -f2 | sort -u > $species/FP/fpSegDupEnsg.txt
join -t$'\t' -1 1 -2 1 -o '2.3 2.4 2.5 1.1' <(cut -f1 $species/TP/tpGeneTranscripts.tsv | sort -u) <(sort $ensemblDir/geneTranscript.ensembl74.tsv) | sort -u | overlapSelect stdin $ucscDir/hg19.segmentalDuplications.bed -selectFmt=bed -statsOutput stdout | tail -n+2 | cut -f2 | sort -u > $species/TP/tpSegDupEnsg.txt

# Some species specific R analysis
if [[ $FPs > 0 && $TPs > 0 ]]
then
echo -e "\n  ...species exons/bases deletions in TP and FP calls (in R)..."
cd $species
sed -e "s/currSpecies/$species/g" $src/R/transcriptsTpFp.R > transcriptsTpFp.R
Rscript transcriptsTpFp.R
cd ../
fi
#sed -e "s/sp_analysis/$species/g" TpFpFeatures.R > $species/TpFpFeatures.R
#Rscript TpFpFeatures.R

# Update the precision-scatter plot (via R) whenever a new species is added
echo -e "\n  ...updating precision by Ensembl74 scatter plot (R)..."
Rscript $src/R/precisionScatter.R
cp callsVsEnsembl.png ../

# Plot TP/FP distributions for various attributes (e.g., fraction of deleted exons...)
#echo -e "\n  ...plotting some distrbutions of various TP/FP attributes (R)..."
#Rscript $src/R/deletionDistributionPlacentals.R 

# creating filter log:
echo -e "$hfiltersApplied\t$ofiltersApplied\t$fracBasesDeleted" > $species/filters.log

echo -e "\n DONE!"
