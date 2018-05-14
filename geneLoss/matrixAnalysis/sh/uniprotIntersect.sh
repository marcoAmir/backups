#!/bin/bash -e

# This script will take the human gene loss calls in a particular species and will check the status of the orthologs in Uniprot

ensemblDir=/cluster/u/amirma/geneLoss/data/Ensembl
uniprotDir=/cluster/u/amirma/geneLoss/data/Uniprot
ncbiDir=/cluster/u/amirma/geneLoss/data/NCBI
speciesList=/cluster/u/amirma/geneLoss/matrixAnalysis/data/speciesList.txt

# Usage
if [ "$#" -ne 2 ]; then
  echo -e "\nUsage:\n$0 geneLossFile species \n"
  exit 1
fi

# Inputs
geneLossCalls=$1
species=$2

assembly=`grep -w $species $speciesList | cut -f1`
species=`echo $species | sed -e "s/_//g"`
taxon=`grep $species $ncbiDir/taxonID.txt | cut -f3`

echo -e "\n\t$species\t$assembly\t$taxon"

geneLosses=`awk -F'\t' '{if($3==1) print $0}' $geneLossCalls | grep -w $assembly | sort -u | wc -l`
echo -e "\t$geneLosses human-genes have been called lost in $assembly\n"

# I - get the ensembl orthologs and infer the intersection with uniprot:
echo -e "I - Uniprot intersection via Ensembl orthologs:"
join -t$'\t' -1 1 -2 1 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10' <(awk -F'\t' '{if($3==1) print $0}' $geneLossCalls | grep -w $assembly | sort -u) <(sort $ensemblDir/OrthoTables74/human${species}.EnsemblOrtho.biomart74.map) > ensTmp
ensOrtho=`cut -f6 ensTmp | sort -u | wc -l`
ensg=`cut -f1 ensTmp | sort -u | wc -l`
protCod=`awk -F'\t' '{IGNORECASE=1; if($9=="protein_coding") print $6}' ensTmp | sort -u | wc -l`
protCodSameName=`awk -F'\t' '{IGNORECASE=1; if($9=="protein_coding" && $3==$8) print $6}' ensTmp | sort -u | wc -l`
echo -e "\t$ensOrtho $species ensembl orthologs were found for $ensg of the $geneLosses human gene losses"
echo -e "\t$protCod out of $ensOrtho $species ensembl orthologs are protein-coding"
echo -e "\t$protCodSameName out of $ensOrtho $species ensembl orthologs are protein-coding with same name as in human"
echo -e "    Inferring uniprot IDs from orthologs..."
join -t$'\t' -1 3 -2 1 -o '1.1 1.2 2.1 2.2 2.3 2.4 2.5 2.7' <(join -t$'\t' -1 2 -2 1 -o '1.1 1.2 2.3' <(cut -f1,6 ensTmp | sort -t$'\t' -k2,2 ) <(sort $ensemblDir/ensemblUniprot/ensemblUniprot.biomart74.${species}.tsv) | awk -F'\t' '{if($3!="") print $0}' | sort -u | sort -t$'\t' -k3,3) <(sort $uniprotDir/uniprot_sprot.parsed.eutheria) | awk -F'\t' '{if($6=='$taxon') print $0}' > uniTmp
ensUni=`cut -f2 uniTmp | sort -u | wc -l`
Uni=`cut -f3 uniTmp | sort -u | wc -l`
echo -e "\t$Uni uniprot identifiers were found for $ensUni of the $ensOrtho $species ensembl orthologs"
echo -e "\tLevel of evidence:"
cut -f8 uniTmp | sed -e "s/;//g;s/://g" | sort | awk '{for(i=2; i<=NF; i++) printf $i" "; print "(tier "$1")"}' | uniq -c | awk '{print "    "$0}'
FPs=`cut -f1 uniTmp | sort -u | wc -l`
echo -e "\tTotal of $FPs human-genes false-positive genes in $geneLosses $species gene-loss calls"
cut -f1 uniTmp | sort -u | awk -F'\t' '{print "\t  "$0}'

# I - get the orthologs from inparanoid  (using uniprot IDs):
echo -e "\nII - Uniprot intersection via inParanoid:"
join -t$'\t' -1 1 -2 1 -o '1.1 2.3' <(awk -F'\t' '{if($3==1) print $0}' $geneLossCalls | grep -w $assembly | sort -u) <(sort $ensemblDir/ensemblUniprot/ensemblUniprot.biomart74.Human.tsv) | awk -F'\t' '{if($2!="") print $0}' | sort -u > ensUniHuman
ensg=`cut -f1 ensUniHuman | sort -u | wc -l`
uni=`cut -f2 ensUniHuman | sort -u | wc -l`
echo -e "\t$uni uniprot identifiers map to $ensg of the $geneLosses ensembl human genes called lost in $species"
join -t$'\t' -1 2 -2 5 -o '1.1 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10' <(sort -t$'\t' -k2 ensUniHuman) <(sort -t$'\t' -k5,5 $uniprotDir/inparanoid/orthoTable.Human-${species}) | sort -u | sort -nk2,2 -nk3,3 > ensUniHuman${species}
uniH=`cut -f6 ensUniHuman${species} | sort -u | wc -l`
uniO=`cut -f10 ensUniHuman${species} | sort -u | wc -l`
clusters=`cut -f2 ensUniHuman${species} | sort -u | wc -l`
clustersPerfect=`awk -F'\t' '{if($5==1 && $7=="100%" && $9==1 && $11=="100%") print $2}' ensUniHuman${species} | sort -u | wc -l`
echo -e "\t$uniH of the $uni uniprot identifiers map by inparanoid to $uniO $species uniprot orthologs in $clusters clusters"
echo -e "\t$clustersPerfect of the $clusters inparanoid orthology clusters perfectly mapped (score=1.0; bootstrap=100%)"
join -t$'\t' -1 10 -2 1 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 2.4 2.7' <(sort -t$'\t' -k10,10 ensUniHuman${species}) <(sort $uniprotDir/uniprot_sprot.parsed.eutheria) | awk -F'\t' '{if($12=='$taxon') print $0}' > ensUniHuman${species}full
uniOevi=`cut -f10 ensUniHuman${species}full | sort -u | wc -l`
echo -e "\t$uniOevi of the $uniO $species uniprot orthologs have explicit annotations in uniprot (i.e., found at ${uniprotDir}/uniprot_sprot.parsed.eutheria)"
echo -e "\tLevel of evidence:"
cut -f13 ensUniHuman${species}full | sed -e "s/;//g;s/://g" | sort | awk '{for(i=2; i<=NF; i++) printf $i" "; print "(tier "$1")"}' | uniq -c | awk '{print "    "$0}'
FPs=`cut -f1 ensUniHuman${species}full | sort -u | wc -l`
echo -e "\tTotal of $FPs human-genes false-positive genes in $geneLosses $species gene-loss calls"
cut -f1 ensUniHuman${species}full | sort -u | awk -F'\t' '{print "\t  "$0}'

FPs_mutual=`comm -1 -2 <(cut -f1 ensUniHuman${species}full | sort -u) <(cut -f1 uniTmp | sort -u) | wc -l`
echo -e "\n$FPs_mutual common false positive calls in the two methods\n\n"

# clean tmp files
rm -rf ensTmp uniTmp ensUniHuman*











