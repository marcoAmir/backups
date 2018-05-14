#!/bin/bash -e 

# I - This script takes an input gene (ensembl gene id or symbol) and a species assembly name and spits out a report of the calls we made for it
# II - it's basically a dumb caller for another script for computing Ka/Ks (/cluster/u/amirma/geneLoss/matrixAnalysis/placentalMammalLosses/src/nonSynSub.sh), and it will just carve out the relevant information for the gene from there
# III - make sure you are logged into dev


ensemblDir=$GENELOSS/data/Ensembl
uniprotDir=$GENELOSS/data/Uniprot
ncbiDir=$GENELOSS/data/NCBI
speciesList=$GENELOSS/data/speciesList.txt
callTable=$GENELOSS/data/hg19.calltable.unique.placentals
src=$GENELOSS/src
KaKsData=/cluster/u/amirma/geneLoss/matrixAnalysis/placentalMammalLosses/KaKs/speciesData

addKaKsPerExon=1

workDir=$PWD

# Usage
if [ "$#" -ne 2 ]; then
  echo -e "\nUsage:\n$0 geneID/symbol species \n"
  exit 1
fi

gene=$1
species=$2

geneSymbol=`grep -w $gene $ensemblDir/humanGeneSymbol.humanEnsembl.biomart74.NoSyn.map | cut -f1`
ensGeneId=`grep -w $gene $ensemblDir/humanGeneSymbol.humanEnsembl.biomart74.NoSyn.map | cut -f2`
assembly=`grep -i -w $species $speciesList | cut -f4`

echo -e "\n$geneSymbol $ensGeneId $species $assembly\n\n\t...Retreiving data...\n"

$src/sh/nonSynSub.sh $assembly

grep -f <(grep -w $ensGeneId $ensemblDir/geneTranscript.ensembl74.tsv | cut -f2 | sort -u) <(cat ${KaKsData}/transcriptNonSynRate.${assembly}.tsv) > tmp

cut -f2 tmp | sort | uniq -c | awk '{if($1>1) print "\nWarning! Detected multiple calls for multiple chains for transcript "$2}'

i=0
while read l
do
	i=$(( i + 1 )) 
	tx=`echo $l | awk '{print $2}'`
	stop=`echo $l | awk '{print $4}'`
	fs=`echo $l | awk '{print $5}'`
	KaKs=`echo $l | awk '{print $11}'`
	grep $assembly $callTable | grep $tx | head -1 > tmp2
	delBases=`awk -F'\t' '{print $10}' tmp2`
	totalBases=`awk -F'\t' '{print $14}' tmp2`
	delExons=`awk -F'\t' '{print $5}' tmp2`
	totalExons=`awk -F'\t' '{print $15}' tmp2`
	if [ $i -eq 1 ]
	then
		echo -e "transcriptID\tstop.codons\tframeshifts\tnum.Deleted.Exons\tnum.Total.Exons\tnum.Deleted.CodingBases\tnum.Total.CodingBases\tKa/Ks" > tmp3
		if [ $addKaKsPerExon -eq 1 ]; then
			echo -e "transcriptID\tstop.codons\tframeshifts\tnum.Deleted.Exons\tnum.Total.Exons\tnum.Deleted.CodingBases\tnum.Total.CodingBases\tKa/Ks\tKa/Ks per exon" > tmp3
		fi
	fi
	if [ $addKaKsPerExon -ne 1 ]; then
		echo -e "$tx\t$stop\t$fs\t$delExons\t$totalExons\t$delBases\t$totalBases\t$KaKs" >> tmp3
	fi
	if [ $addKaKsPerExon -eq 1 ]; then
		KaKsPerExon=`grep $tx $callTable | grep $assembly | cut -f17 | sed -e "s/,/\n/g" | awk -F':' '{split($0,a,":"); if(a[2]!=0 && a[4]!=0 && a[3]/a[4]!=0){if(NR==1) {printf "{%.2f", (a[1]/a[2])/(a[3]/a[4])}; if(NR>1) {printf "; %.2f", (a[1]/a[2])/(a[3]/a[4])}} else {if(NR==1) {printf "{NA"}; if(NR>1) {printf "; NA"}}} END {print "}"}'`
		echo -e "$tx\t$stop\t$fs\t$delExons\t$totalExons\t$delBases\t$totalBases\t$KaKs\t$KaKsPerExon" >> tmp3
	fi
done < tmp

cat tmp3 | column -t -s$'\t'

#awk -F'\t' '{if(NR==1) {print "transcriptID\tstop.codons\tframeshifts\tnum.Deleted.Exons\tnum.Total.Exons\tnum.Deleted.CodingBases\tnum.Total.CodingBases\tKa/Ks (transcript)\ttranscript LoF call"}; print $2"   \t"$4"          \t"$5"        \t"$6"                 \t"$7"        \t"$8"      \t"$9"      \t"$11"                       \t"$10}' tmp | column -t

rm -rf tmp*

echo ""

#awk -F'\t' '{printf "\t\t"$2}' tmp | awk '{print $0}'
#awk -F'\t' '{if(NR==1 ) {printf "stop codons:      "} else {printf "\t"$4}}' tmp | awk '{print $0}'
#awk -F'\t' '{if(NR==1 ) {printf "frameshifts:      "} else {printf "\t"$5}}' tmp | awk '{print $0}'
#awk -F'\t' '{if(NR==1 ) {printf "num.Deleted.Exons:"} else {printf "\t"$6}}' tmp | awk '{print $0}'
#awk -F'\t' '{if(NR==1 ) {printf "num.Total.Exons:"} else {printf "\t"$7}}' tmp | awk '{print $0}'
#awk -F'\t' '{if(NR==1 ) {printf "Ka/Ks (transcript):"} else {printf "\t"$11}}' tmp | awk '{print $0}'
#awk -F'\t' '{if(NR==1 ) {printf "transcript LoF call:"} else {printf "\t"$8}}' tmp | awk '{print $0}'



