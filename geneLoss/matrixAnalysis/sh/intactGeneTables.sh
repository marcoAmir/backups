#!/bin/bash -e

# Here, for a given species, I'll identify the intact protein-coding genes that are stable in ensembl 
ensembl=$GENELOSS/data/Ensembl
ucsc=$GENELOSS/data/UCSC
curatedGenes=$GENELOSS/data/curatedGeneLosses.validated.txt

# Usage
if [ "$#" -ne 2 ]; then
  echo -e "\nUsage:\n$0 species annotatedByEnsembl(1/0)\n"
  exit 1
fi

species=$1
assembly=`grep $species $GENELOSS/data/speciesList.txt | cut -f4`
if [ "$species" = "Cat" ]
then
assembly=felCat5
fi
species=`echo $species | sed -e "s/_//g"`
logFile=${species}Human.log

useCanonicals=1
filter_SegDups_Deprecated=1	# if selected I'll remove all genes that intersect segmental-duplications & deprecated genes (genes that are protein-coding in Ensembl74 but not in Ensembl78)
annotatedByEnsembl=$2		# if the species is NOT annotated by ensembl, I'll just collect the features and partition to canonical and non-canonical human transcripts

if [ $annotatedByEnsembl -ne 0 ]
then
## I'll only take those with evidence from the literature
#awk -F'\t' '{if($4~/[0-9]/) print $0}' $curatedGenes | grep $species | cut -f1 | sort -u > ${species}.curated
awk -F'\t' '{if($1~/ENSG/) print $0}' $curatedGenes | grep $species | cut -f1 | sort -u > ${species}.curated

echo -e "\n"$species"\t"$assembly
echo -e "\n"$species"\t"$assembly > $logFile

# Step 1 - get the transcripts of the given species that maintain 'protein_coding' biotype across multiple ensembl version (Specifically, from ens74 to ens81):
ens74trans=`awk -F'\t' '{if($6=="protein_coding") print $2}' $ensembl/speciesTranscripts/geneTranscript.${species}.ensembl74.tsv | sort -u | wc -l`
ens81trans=`awk -F'\t' '{if($6=="protein_coding") print $2}' $ensembl/speciesTranscripts/geneTranscript.${species}.ensembl81.tsv | sort -u | wc -l`
join -t$'\t' -1 2 -2 2 -o '1.1 1.2 1.3 1.4 1.5 1.6 2.1 2.2 2.3 2.4 2.5 2.6' <(sort -k2,2 $ensembl/speciesTranscripts/geneTranscript.${species}.ensembl74.tsv) <(sort -k2,2 $ensembl/speciesTranscripts/geneTranscript.${species}.ensembl81.tsv) | awk -F'\t' '{if($6=="protein_coding" && $12=="protein_coding" && $3==$9) print $0}' | cut -f1,2 | sort -u > ${species}.intact.proteinCodingTx
ensOverlap=`cut -f2 ${species}.intact.proteinCodingTx | sort -u | wc -l`

# Step 2 - intersect the list from step 1 with a list of complete transcript model in mouse
if [ "$species" = "Dolphin" ]	# pretty annoying but i don't have the ensGene table for turTru2
then
assembly=turTru1
fi
if [ "$species" = "Alpaca" ]   # pretty annoying but i don't have the ensGene table for vicPac2
then
assembly=vicPac1
fi
join -t$'\t' -1 2 -2 1 -o '1.1 1.2' <(sort -k2 ${species}.intact.proteinCodingTx) <(hgsql ${assembly} -Ne "SELECT name, chrom, txStart, txEnd FROM ensGene WHERE cdsStartStat='cmpl' AND cdsEndStat='cmpl'" | sort) | sort -u > tmp; mv tmp ${species}.intact.proteinCodingTx
if [ "$species" = "Dolphin" ]	# after the sql query I'll just revert the assembly back to turTru2
then
assembly=turTru2
fi
if [ "$species" = "Alpaca" ]   # pretty annoying but i don't have the ensGene table for vicPac2
then
assembly=vicPac2
fi
ensIntacts=`wc -l ${species}.intact.proteinCodingTx | cut -d" " -f1`
echo -e "\n\t$ens74trans   $species protein-coding transcripts in Ensembl biomart 74\n\t$ens81trans   $species protein-coding transcripts in Ensembl biomart 81\n\t$ensOverlap   $species shared protein-coding transcripts\n\t$ensIntacts   $species shared protein-coding (ens74 to ens81) and complete models" >> $logFile
echo -e "\n\t$ens74trans   $species protein-coding transcripts in Ensembl biomart 74\n\t$ens81trans   $species protein-coding transcripts in Ensembl biomart 81\n\t$ensOverlap   $species shared protein-coding transcripts\n\t$ensIntacts   $species shared protein-coding (ens74 to ens81) and complete models" 

# Step 3 - Now I'll expand the list of intact transcript (gene, transcript pairs) that passed the filters in steps 1,2 with human orthologs at the transcript level.
# Namely the transcript IDs of the query species (e.g., mouse ENSMUST do not play a role from here on)
# I will use the Ensembl74 tables
# Output from this step is a table mapping species gene identifiers (intact protein-coding genes by steps 1,2) to human gene idenftiers along with a list of their protein-coding transcripts
join -t$'\t' -1 1 -2 6 -o '1.1 2.1 2.2 2.3 2.4' <(sort ${species}.intact.proteinCodingTx) <(sort -k6,6 $ensembl/OrthoTables74/human${species}.EnsemblOrtho.biomart74.map) | sort -u | awk -F'\t' '{if($5=="protein_coding") print $1"\t"$2"\t"$3"\t"$4}' > ${species}.intact.proteinCodingTx.humanOrtho
# I will remove here ambigous gene-gene orthologs where a query species gene is mapped to more than a single human gene, or vice versa where a human-gene is mapped to more than a single gene in the query species
cut -f1,2 ${species}.intact.proteinCodingTx.humanOrtho | sort -u | cut -f2 | sort | uniq -c | awk '{if($1>1) print $2}' | sort -u > nonUniqeHumanGenes.txt
cut -f1,2 ${species}.intact.proteinCodingTx.humanOrtho | sort -u | cut -f1 | sort | uniq -c | awk '{if($1>1) print $2}' | sort -u > nonUniqe${species}Genes.txt
grep -v -f nonUniqe${species}Genes.txt <(grep -v -f nonUniqeHumanGenes.txt ${species}.intact.proteinCodingTx.humanOrtho) > tmp; mv tmp ${species}.intact.proteinCodingTx.humanOrtho 


# Step 4 - remove incomplete human transcripts
join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4' <(hgsql hg19 -Ne "SELECT name FROM ensGene WHERE cdsStartStat='cmpl' AND cdsEndStat='cmpl'" | sort) <(sort -k3 ${species}.intact.proteinCodingTx.humanOrtho) > tmp ; mv tmp ${species}.intact.proteinCodingTx.humanOrtho

# Step 5 - Now, expand with:
		# chains (2.3), stop-codons (2.4), frameshifts (2.11), number of deleted exons (2.5), number of total exons (2.15), number of deleted bases (2.10), number of deleted bases (2.14), percent.id (2.16)
join -t$'\t' -1 3 -2 2 -o '2.1 1.1 1.2 1.3 1.4 2.3 2.4 2.11 2.5 2.15 2.10 2.14 2.16' <(sort -t$'\t' -k3,3 ${species}.intact.proteinCodingTx.humanOrtho) <(grep ${assembly} $GENELOSS/data/hg19.calltable.unique.placentals | sort -t$'\t' -k2,2) | sort -t$'\t' -k4,4 > tmp

# And add the Ka/Ks ratio (pre-computed at: $GENELOSS/data/KaKs)
join -t$'\t' -1 4 -2 2 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 2.2 2.3 2.11' <(sort -t$'\t' -k4,4 tmp) <(sort -t$'\t' -k2,2 $GENELOSS/data/KaKs/speciesData/transcriptNonSynRate.${assembly}.tsv) | awk -F'\t' '{if($6==$15) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9"\t"$10"\t"$9/$10"\t"$11"\t"$12"\t"$11/$12"\t"$13"\t"$16}' | sort -u | sort -k2,2 -k3,3 -k4,4 > ${species}.intact.proteinCodingTx.humanOrtho

# If selected, exclude all the non canonical transcripts
if [ $useCanonicals -ne 0 ]
then
echo -e "\n  Keeping only canonical human transcripts...\n"
join -t$'\t' -1 4 -2 1 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15' <(sort -t$'\t' -k4,4 ${species}.intact.proteinCodingTx.humanOrtho) <(sort $ensembl/hg19.canonicalTranscripts.biomart74) | sort -u | sort -k2,2 -k3,3 -k4,4 > ${species}.intact.proteinCodingTx.canonical.humanOrtho
fi


rm -rf tmp

ensHuman=`cut -f4 ${species}.intact.proteinCodingTx.humanOrtho | sort -u | wc -l`
ensg=`cut -f3 ${species}.intact.proteinCodingTx.humanOrtho | sort -u | wc -l`
qgenes=`cut -f2 ${species}.intact.proteinCodingTx.humanOrtho | sort -u | wc -l`
echo -e "Summary\n=======\n\t$ensHuman   human transcripts were found in $ensg genes orthologues to $qgenes intact genes in ${species}" >> $logFile
echo -e "Summary\n=======\n\t$ensHuman   human transcripts were found in $ensg genes orthologues to $qgenes intact genes in ${species}" 

# Step 6 - now take all the human transcript that do not belong to genes with orthology to the intact genes in the query species
join -t$'\t' -1 2 -2 2 -o '1.1 2.1 2.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10' <(join -t$'\t' -1 2 -2 2 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 2.2 2.3 2.11' <(join -t$'\t' -1 1 -2 2 -o '2.1 2.2 2.3 2.4 2.11 2.5 2.15 2.10 2.14 2.16' -v 2 <(cut -f4 ${species}.intact.proteinCodingTx.humanOrtho | sort -u) <(grep $assembly $GENELOSS/data/hg19.calltable.unique.placentals | sort -t$'\t' -k2,2) | sort -t$'\t' -k2,2) <(sort -t$'\t' -k2,2 $GENELOSS/data/KaKs/speciesData/transcriptNonSynRate.${assembly}.tsv) | awk -F'\t' '{if($2==$11 && $3==$12) print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$13}' | sort -t$'\t' -k2,2) <(sort -t$'\t' -k2,2 $ensembl/geneTranscript.ensembl74.tsv) | sort -u > tmp
join -t$'\t' -1 2 -2 2 -o '1.1 1.2 1.3 2.1 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11' <(sort -t$'\t' -k2,2 tmp) <(sort -t$'\t' -k2,2 $ensembl/humanGeneSymbol.humanEnsembl.biomart74.NoSyn.map) | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$7/$8"\t"$9"\t"$10"\t"$9/$10"\t"$11"\t"$12}' | sort -u > ${species}.remainingCalls.proteinCodingTx.humanOrtho
join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10 2.11 2.12 2.13 2.14' <(hgsql hg19 -Ne "SELECT name FROM ensGene WHERE cdsStartStat='cmpl' AND cdsEndStat='cmpl'" | sort) <(sort -t$'\t' -k3,3 ${species}.remainingCalls.proteinCodingTx.humanOrtho) | sort -u > tmp; mv tmp ${species}.remainingCalls.proteinCodingTx.humanOrtho

# If selected, exclude all the non canonical transcripts
if [ $useCanonicals -ne 0 ]
then
echo -e "\n  Keeping only canonical human transcripts...\n"
join -t$'\t' -1 3 -2 1 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14' <(sort -t$'\t' -k3,3 ${species}.remainingCalls.proteinCodingTx.humanOrtho) <(sort $ensembl/hg19.canonicalTranscripts.biomart74) | sort -u | sort -k2,2 -k3,3 -k4,4 > ${species}.remainingCalls.proteinCodingTx.canonical.humanOrtho
fi

RemEnsHuman=`cut -f2 ${species}.remainingCalls.proteinCodingTx.humanOrtho | sort -u | wc -l`
RemEnstHuman=`cut -f3 ${species}.remainingCalls.proteinCodingTx.humanOrtho | sort -u | wc -l`
echo -e "\t$RemEnstHuman   human transcripts were found in $RemEnsHuman genes unmapped to intact protein-coding orthologues genes in ${species}" >> $logFile
echo -e "\t$RemEnstHuman   human transcripts were found in $RemEnsHuman genes unmapped to intact protein-coding orthologues genes in ${species}"

# Step 7 - Same as in 6 but here we will select only the trancript which in we called lost by the selected threshold
join -t$'\t' -1 2 -2 2 -o '1.1 2.1 2.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10' <(join -t$'\t' -1 2 -2 2 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 2.2 2.3 2.11' <(join -t$'\t' -1 1 -2 2 -o '2.1 2.2 2.3 2.4 2.11 2.5 2.15 2.10 2.14 2.16' -v 2 <(cut -f4 ${species}.intact.proteinCodingTx.humanOrtho | sort -u) <(grep $assembly $GENELOSS/data/hg19.calltable.unique.placentals | sort -t$'\t' -k2,2) | sort -t$'\t' -k2,2) <(sort -t$'\t' -k2,2 $GENELOSS/data/KaKs/speciesData/transcriptNonSynRate.${assembly}.tsv) | awk -F'\t' '{if($2==$11 && $3==$12) print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$13}' | sort -t$'\t' -k2,2) <(sort -t$'\t' -k2,2 $ensembl/geneTranscript.ensembl74.tsv) | sort -u > tmp
join -t$'\t' -1 2 -2 2 -o '1.1 1.2 1.3 2.1 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11' <(sort -t$'\t' -k2,2 tmp) <(sort -t$'\t' -k2,2 $ensembl/humanGeneSymbol.humanEnsembl.biomart74.NoSyn.map) | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$7/$8"\t"$9"\t"$10"\t"$9/$10"\t"$11"\t"$12}' | sort -u > ${species}.lossCalls.proteinCodingTx.humanOrtho
join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10 2.11 2.12 2.13 2.14' <(hgsql hg19 -Ne "SELECT name FROM ensGene WHERE cdsStartStat='cmpl' AND cdsEndStat='cmpl'" | sort) <(sort -t$'\t' -k3,3 ${species}.lossCalls.proteinCodingTx.humanOrtho) | sort -u > tmp; mv tmp ${species}.lossCalls.proteinCodingTx.humanOrtho

# If selected, exclude all the non canonical transcripts
if [ $useCanonicals -ne 0 ]
then
#echo -e "\n  Keeping only canonical human transcripts...\n"
join -t$'\t' -1 3 -2 1 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14' <(sort -t$'\t' -k3,3 ${species}.lossCalls.proteinCodingTx.humanOrtho) <(sort $ensembl/hg19.canonicalTranscripts.biomart74) | sort -u | sort -k2,2 -k3,3 -k4,4 > ${species}.lossCalls.proteinCodingTx.canonical.humanOrtho
fi

lostEnsHuman=`cut -f2 ${species}.lossCalls.proteinCodingTx.humanOrtho | sort -u | wc -l`
lostEnstHuman=`cut -f3 ${species}.lossCalls.proteinCodingTx.humanOrtho | sort -u | wc -l`
echo -e "\t$lostEnstHuman   human transcripts were found in $lostEnsHuman genes with Loss-of-Function calls in ${species}\n" >> $logFile
echo -e "\t$lostEnstHuman   human transcripts were found in $lostEnsHuman genes with Loss-of-Function calls in ${species}\n"

# Step 8 - if 'filter_SegDups_Deprecated' is selected create the subset files as well
if [ $filter_SegDups_Deprecated -ne 0 ]
then
cp $ensembl/hg19.transcripts.bed a.bed
cp $ucsc/hg19.segmentalDuplications.bed b.bed
join -t$'\t' -1 1 -2 1 -o '2.2' <(join -t$'\t' -1 1 -2 2 -o '2.1' <(overlapSelect a.bed b.bed -statsOutput stdout | awk -F'\t' '{if($4>0) print $0}' | cut -f2 | sort -u) <(sort -t$'\t' -k2,2 $ensembl/geneTranscript.ensembl74.tsv) | sort -u) <(sort -t$'\t' $ensembl/geneTranscript.ensembl74.tsv) | sort -u > transcriptsInSegDup.txt
# now clean up the files from seg-dups and deprecated genes:
join -t$'\t' -1 1 -2 3 -v 2 <(sort -u $ensembl/genes74noLongerProteinCoding.txt) <(join -t$'\t' -1 1 -2 4 -v 2 <(sort transcriptsInSegDup.txt) <(sort -t$'\t' -k4,4 ${species}.intact.proteinCodingTx.humanOrtho) | sort -u | sort -t$'\t' -k3,3) | sort -u | sort -t$'\t' -k3,3 > ${species}.intact.proteinCodingTx.humanOrtho.filter_SegDups_Deprecated
join -t$'\t' -1 1 -2 3 -v 2 <(sort -u $ensembl/genes74noLongerProteinCoding.txt) <(join -t$'\t' -1 1 -2 4 -v 2 <(sort transcriptsInSegDup.txt) <(sort -t$'\t' -k4,4 ${species}.intact.proteinCodingTx.canonical.humanOrtho) | sort -u | sort -t$'\t' -k3,3) | sort -u | sort -t$'\t' -k3,3 > ${species}.intact.proteinCodingTx.canonical.humanOrtho.filter_SegDups_Deprecated
join -t$'\t' -1 1 -2 3 -v 2 <(sort -u $ensembl/genes74noLongerProteinCoding.txt) <(join -t$'\t' -1 1 -2 4 -v 2 <(sort transcriptsInSegDup.txt) <(sort -t$'\t' -k4,4 ${species}.remainingCalls.proteinCodingTx.humanOrtho) | sort -u | sort -t$'\t' -k3,3) | sort -u | sort -t$'\t' -k3,3 > ${species}.remainingCalls.proteinCodingTx.humanOrtho.filter_SegDups_Deprecated
join -t$'\t' -1 1 -2 3 -v 2 <(sort -u $ensembl/genes74noLongerProteinCoding.txt) <(join -t$'\t' -1 1 -2 4 -v 2 <(sort transcriptsInSegDup.txt) <(sort -t$'\t' -k4,4 ${species}.remainingCalls.proteinCodingTx.canonical.humanOrtho) | sort -u | sort -t$'\t' -k3,3) | sort -u | sort -t$'\t' -k3,3 > ${species}.remainingCalls.proteinCodingTx.canonical.humanOrtho.filter_SegDups_Deprecated
rm -rf a.bed b.bed transcriptsInSegDup.txt
fi

# Step 9 - generate plots in R: (note that these plots show everything including seg-dups transcripts)
sed -e "s/QUERY/\'${species}\'/g" $GENELOSS/src/R/intactGeneDistributions.R > intactGeneDistributions.R
sed -e "s/QUERY/\'${species}\'/g" $GENELOSS/src/R/intactGeneDistributions.clusterAnalysis.R > intactGeneDistributions.clusterAnalysis.R
if [ $useCanonicals -ne 0 ]
then
sed -e "s/QUERY/\'${species}\'/g;s/proteinCodingTx.humanOrtho/proteinCodingTx.canonical.humanOrtho/g;s/\.png/.canonical.png/g;s/\.pdf/.canonical.pdf/g" $GENELOSS/src/R/intactGeneDistributions.R > intactGeneDistributions.R
sed -e "s/QUERY/\'${species}\'/g;s/proteinCodingTx.humanOrtho/proteinCodingTx.canonical.humanOrtho/g;s/\.png/.canonical.png/g;s/\.pdf/.canonical.pdf/g;s/\.Quantiles/.canonical.Quantiles/g" $GENELOSS/src/R/intactGeneDistributions.clusterAnalysis.R > intactGeneDistributions.clusterAnalysis.R
fi
chmod 755 intactGeneDistributions.R
Rscript intactGeneDistributions.clusterAnalysis.R
Rscript intactGeneDistributions.R
fi

if [ $annotatedByEnsembl -eq 0 ]
then
grep ${assembly} $GENELOSS/data/hg19.calltable.placentals | awk -F'\t' '{print $1"\t"$2"\t"$4"\t"$11"\t"$5"\t"$15"\t"$5/$15"\t"$10"\t"$14"\t"$10/$14"\t"$16}' > tmp1
join -t$'\t' -1 2 -2 2 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 2.11' <(sort -t$'\t' -k2,2 tmp1) <(sort -t$'\t' -k2,2 $GENELOSS/data/KaKs/speciesData/transcriptNonSynRate.${assembly}.tsv) > tmp2
join -t$'\t' -1 2 -2 2 -o '1.1 1.2 1.3 2.1 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13' <(join -t$'\t' -1 2 -2 2 -o '1.1 2.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12' <(sort -t$'\t' -k2,2 tmp2) <(sort -t$'\t' -k2,2 $ensembl/geneTranscript.ensembl74.tsv) | sort -t$'\t' -k2,2) <(sort -t$'\t' -k2,2 $ensembl/humanGeneSymbol.humanEnsembl.biomart74.NoSyn.map) | sort -u | sort -t$'\t' -k2,2 | awk -F'\t' '{print $1"\tNA\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' > tmp3
if [ $filter_SegDups_Deprecated -ne 0 ]
then
cp $ensembl/hg19.transcripts.bed a.bed
cp $ucsc/hg19.segmentalDuplications.bed b.bed
join -t$'\t' -1 1 -2 1 -o '2.2' <(join -t$'\t' -1 1 -2 2 -o '2.1' <(overlapSelect a.bed b.bed -statsOutput stdout | awk -F'\t' '{if($4>0) print $0}' | cut -f2 | sort -u) <(sort -t$'\t' -k2,2 $ensembl/geneTranscript.ensembl74.tsv) | sort -u) <(sort -t$'\t' $ensembl/geneTranscript.ensembl74.tsv) | sort -u > transcriptsInSegDup.txt
join -t$'\t' -1 1 -2 3 -v 2 <(sort -u $ensembl/genes74noLongerProteinCoding.txt) <(join -t$'\t' -1 1 -2 4 -v 2 <(sort transcriptsInSegDup.txt) <(sort -t$'\t' -k4,4 tmp3) | sort -u | sort -t$'\t' -k3,3) | sort -u | sort -t$'\t' -k3,3 > ${species}.all.proteinCodingTx.humanOrtho.filter_SegDups_Deprecated
cp ${species}.all.proteinCodingTx.humanOrtho.filter_SegDups_Deprecated tmp3
fi
if [ $useCanonicals -ne 0 ]
then
join -t$'\t' -1 4 -2 1 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15' <(sort -t$'\t' -k4,4 tmp3) <(sort $ensembl/hg19.canonicalTranscripts.biomart74) | sort -u | sort -t$'\t' -k3,3 > ${species}.all.proteinCodingTx.canonical.humanOrtho
if [ $filter_SegDups_Deprecated -ne 0 ]
then 
mv ${species}.all.proteinCodingTx.canonical.humanOrtho ${species}.all.proteinCodingTx.canonical.humanOrtho.filter_SegDups_Deprecated
fi
fi
fi
