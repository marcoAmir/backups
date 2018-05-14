#!/bin/bash -e

ensemblDir=/cluster/u/amirma/geneLoss/data/Ensembl
callTable=/cluster/u/amirma/geneLoss/matrixAnalysis/data/hg19.calltable
pairwiseBranchLengths=/cluster/u/amirma/geneLoss/data/UCSC/multiz100way/pairWiseDistances.tsv
KaKsData=/cluster/u/amirma/geneLoss/matrixAnalysis/placentalMammalLosses/KaKs/speciesData

if [ "$#" -ne 1 ]; then
  echo -e "\nUsage:\n  $0 gene (symbol or id)\n"
  exit 1
fi

gene=$1

geneSymbol=`grep -w $gene $ensemblDir/humanGeneSymbol.humanEnsembl.biomart74.NoSyn.map | cut -f1`
geneID=`grep -w $gene $ensemblDir/humanGeneSymbol.humanEnsembl.biomart74.NoSyn.map | cut -f2`
echo -e "\n\t$geneSymbol\t$geneID\n"

grep $geneID $ensemblDir/geneTranscript.ensembl74.tsv | cut -f2 | sort -u > currTranscripts

for f in $KaKsData/transcript* 
do 
	echo "processing $f"
	join -t$'\t' -1 1 -2 2 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10 2.11 2.12 2.13' <(sort currTranscripts) <(sort -t$'\t' -k2,2 $f) | sort -u > tmp
	nChains=`cut -f3 tmp | sort -u | wc -l`
	if [[ $nChains -gt 1 ]]
	then
		echo -e "\n$f ignored: ambigious pairwise alignement with $nChains chains - resolve manually\n"
		echo -e "$f" >> ${geneSymbol}_unresolvedPairwise
	else
		cat tmp >> $geneSymbol
	fi
done
