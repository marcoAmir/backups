#!/bin/bash -e

intactGeneDir=/cluster/u/amirma/geneLoss/intactGenes

for species in Cat Chimp Dog Dolphin Elephant Gorilla GuineaPig Horse Marmoset Megabat Microbat Mouse Panda Rabbit Rat
do
echo -e $species
# canonical transcripts
cut -f1,3-15 $intactGeneDir/${species}.intact.proteinCodingTx.canonical.humanOrtho | awk -F'\t' '{print $0"\tintact\t0"}' > ${species}.allFeatures
cat $intactGeneDir/${species}.remainingCalls.proteinCodingTx.canonical.humanOrtho | awk -F'\t' '{print $0"\tother\t1"}'>> ${species}.allFeatures
join -t$'\t' -1 2 -2 1 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 2.1' -a 1 <(sort -t$'\t' -k2,2 ${species}.allFeatures) <(sort $intactGeneDir/${species}.curated) | awk -F'\t' '{if($17~/ENSG/) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\tcurated\t2"} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}}' > tmp; mv tmp ${species}.allFeatures
# all the transcript
cut -f1,3-15 $intactGeneDir/${species}.intact.proteinCodingTx.humanOrtho | awk -F'\t' '{print $0"\tintact\t0"}' > ${species}.allTx.allFeatures
cat $intactGeneDir/${species}.remainingCalls.proteinCodingTx.humanOrtho | awk -F'\t' '{print $0"\tother\t1"}'>> ${species}.allTx.allFeatures
join -t$'\t' -1 2 -2 1 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 2.1' -a 1 <(sort -t$'\t' -k2,2 ${species}.allTx.allFeatures) <(sort $intactGeneDir/${species}.curated) | awk -F'\t' '{if($17~/ENSG/) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\tcurated\t2"} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}}' > tmp; mv tmp ${species}.allTx.allFeatures
done
