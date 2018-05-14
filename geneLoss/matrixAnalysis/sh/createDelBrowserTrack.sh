#!/bin/bash -e

echo -e '\n\tmake sure you are in dev\n'

assembly=$1
species=`grep $assembly /cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/data/speciesList.txt | cut -f2`

cut -f1-3 ${assembly}.deletions.bed > current.bed
hgLoadBed hg19 bejamDeletions${assembly} current.bed -allowStartEqualEnd

echo -e "\n#### COPY THE FOLLOWING BLOCK INTO AN APPROPRIATE PLACE IN /cluster/u/amirma/git/browserTracks/users/amirma/data/human/hg19/trackDb.ra \n"
echo -e "\ntrack bejamDeletions${assembly}"
echo -e "shortLabel Bej ${assembly} LargeDel"
echo -e "longLabel ${species} Predicted Deletions in Coding Regions> 20 Bases"
echo -e "priority 1"
echo -e "group compGeno"
echo -e "visibility hide"
echo -e "color 100,100,100"
echo -e "type bed 3 ."




