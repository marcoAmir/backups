#!/bin/bash -e

src=/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/src/
data=/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data/

excludeSegDup=1

rm -rf tmp*

if [ -f "refTranscriptsGene.txt" ]
then
	echo "refTranscriptsGene.txt found."
else
	join -t$'\t' -1 4 -2 2 -o '1.4 2.1' <(sort -t$'\t' -k4,4 refTranscripts.bed) <(sort -t$'\t' -k2,2 /cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data/Ensembl/geneTranscript.ensembl74.tsv) | sort -u > refTranscriptsGene.txt
fi

for assembly in ailMel1 bosTau7 calJac3 camFer1 canFam3 capHir1 cavPor3 cerSim1 chiLan1
do
	echo $assembly
	join -t$'\t' -1 2 -2 1 -o '1.1 2.2 2.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8' <(overlapSelect refTranscripts.bed ${assembly}.deletions.bed -statsOutput stdout | sort -t$'\t' -k2,2) <(sort refTranscriptsGene.txt) | cut -f1-3 | sort -u  > tmpDelGeneTx
	while read l
	do 
		delID=`echo $l | cut -d" " -f4`
		#echo $delID
		grep $delID tmpDelGeneTx | cut -f2 | sort -u > tmpGene
		while read currGene
		do
			nMissedTxs=`comm -2 -3 <(grep $currGene refTranscriptsGene.txt | cut -f1 | sort -u) <(grep $delID tmpDelGeneTx | cut -f3 | sort -u) | wc -l`
			if [ $nMissedTxs -eq 0 ]
			then
				grep $delID ${assembly}.deletions.bed >> ${assembly}.deletions.allTx.bed 
			fi
		done < tmpGene 
	done < ${assembly}.deletions.bed
	sort -u ${assembly}.deletions.allTx.bed > tmp; mv tmp ${assembly}.deletions.allTx.bed
	ex=`comm -1 -3 <(sort -u ${assembly}.deletions.bed) <(sort -u ${assembly}.deletions.allTx.bed) | wc -l`
	if [ $ex -eq 0 ]
	then
		mv ${assembly}.deletions.allTx.bed ${assembly}.deletions.bed
	else
		echo -e "\n...Warning with $assembly:\n\t\t${assembly}.deletions.allTx.bed not a subset of ${assembly}.deletions.bed \n"
	fi
	if [ $excludeSegDup -eq 1 ]
	then
		echo -e "...Excluding deletions that overlap segmental duplications..."
		overlapSelect $data/UCSC/hg19.segmentalDuplications.bed ${assembly}.deletions.bed -nonOverlapping stdout > tmp; mv tmp ${assembly}.deletions.bed
	fi
done
