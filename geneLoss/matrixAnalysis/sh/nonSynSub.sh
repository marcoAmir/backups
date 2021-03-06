#!/bin/bash -e

# make sure you are in dev

ensemblDir=$GENELOSS/data/Ensembl
callTable=$GENELOSS/data/hg19.calltable.unique.placentals
pairwiseBranchLengths=$GENELOSS/data/UCSC/multiz100way/pairWiseDistances.tsv
KaKsDir=$GENELOSS/data/KaKs/speciesData

if [ "$#" -ne 1 ]; then
  echo -e "\nUsage:\n  $0 assembly\n"
  exit 1
fi

assembly=$1
a=`echo $assembly | sed -e "s/[0-9]//g"`

#geneCalls=$2

#exonDelThreshold=$3

bl=`grep -w $a $pairwiseBranchLengths | cut -f3`		# branch length from human

# the following command will write:
# 1- species; 2- transcriptID; 3- chain; 4 - premature stop; 5 - frameshifts; 6  -num.deleted.exons; 7 - num total transcript exons; 8 - call (1-loss/ 0 /?) ; 9 - Ka non-syn subst (mean over transcript exons) ; 10 - Ks - syn subst. ; 11 - Ka/Ks ratio ; 12 - non-synonymous substituion rate normalized by the number of transcript exons;  13- synonymous substituion rate normalized by the branch length

echo -e "assembly\tenstID\tchainID\tstop\tframeshifts\tdelExons\ttotalExons\tdelCodingBases\ttotalCodingBases\tcall\tKa\tKs\tKaKsRatio\tnormedKa\tnormedKs" > $KaKsDir/transcriptNonSynRate.${assembly}.tsv

## For only complete transcripts:
#join -t$'\t' -1 2 -2 1 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11' <(awk -F'\t' '{if($1~/'$assembly'/) print $1"\t"$2"\t"$3"\t"$4"\t"$11"\t"$5"\t"$15"\t"$10"\t"$14"\t"$24"\t"$17}' $callTable | sort -t$'\t' -k2,2) <(hgsql hg19 -Ne "SELECT name FROM ensGene WHERE cdsStartStat='cmpl' AND cdsEndStat='cmpl'" | sort) | tr ',' '\t' | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10; ka=0; ks=0; for(i=11; i<=NF; i++){split($i,a,":"); if(a[2]!=0 && a[1]<a[2]) {b=a[1]/a[2]; ka+=b}; if(a[4]!=0 && a[3]<a[4]) {b=a[3]/a[4]; ks+=b}}; print "\t"ka/$7"\t"ks/$7"\t"(ka/$7)/'$bl'"\t"(ks/$7)/'$bl' }' | awk -F'\t' '{if($12!=0) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$11/$12"\t"$13"\t"$14} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\tNA\t"$13"\t"$14}}' | sort -u >> $KaKsDir/transcriptNonSynRate.${assembly}.tsv

# for all transcripts:
awk -F'\t' '{if($1~/'$assembly'/) print $1"\t"$2"\t"$3"\t"$4"\t"$11"\t"$5"\t"$15"\t"$10"\t"$14"\t"$24"\t"$17}' $callTable | tr ',' '\t' | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10; ka=0; ks=0; for(i=11; i<=NF; i++){split($i,a,":"); if(a[2]!=0 && a[1]<a[2]) {b=a[1]/a[2]; ka+=b}; if(a[4]!=0 && a[3]<a[4]) {b=a[3]/a[4]; ks+=b}}; print "\t"ka/$7"\t"ks/$7"\t"(ka/$7)/'$bl'"\t"(ks/$7)/'$bl' }' | awk -F'\t' '{if($12!=0) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$11/$12"\t"$13"\t"$14} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\tNA\t"$13"\t"$14}}' | sort -u >> $KaKsDir/transcriptNonSynRate.${assembly}.tsv


#awk -F'\t' '{if($1~/'$assembly'/ && $5>'${exonDelThreshold}') print $0}' $callTable > transcriptWithExonDel.${assembly}.tsv
#awk -F'\t' '{if($1~/'$assembly'/ && $5<='${exonDelThreshold}') print $0}' $callTable > transcriptNoExonDel.${assembly}.tsv












