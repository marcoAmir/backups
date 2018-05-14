#!/bin/bash -e 

data=$GENELOSS/data

assembly=mm10

zcat $data/GENCODE/gencode.vM4.2wayconspseudos.gtf.gz | awk -F'\t' '{if(NF>5) print $9}' | awk -F';' '{print $13"\t"$2"|"$3"|"$6"|"$9"|"$10}' | sed -e "s/ parent_id \"//g;s/\"\t/\t/g" | sort > mousePseudogenes.txt

nPsgs=`cut -f1 mousePseudogenes.txt | sort -u | wc -l`

echo -e "\n\t$nPsgs mouse psedogenes found in GENCODE"


# Step 1 - map the mouse pseudogenes identifiers (ENSMUSG) to human gene identifiers (ENSG) via the ensembl orthology tables
# First, take the pairs that are in agreement between differetn ensembl versions (Ens74 to Ens80)
# Then, exclude ambigious mappings
comm -1 -2 <(join -t$'\t' -1 1 -2 6 -o '1.1 2.1' <(cut -f1 mousePseudogenes.txt | sort -u) <(sort -t$'\t' -k6,6 $data/Ensembl/OrthoTables74/humanMouse.EnsemblOrtho.biomart74.map) | sort -u) <(join -t$'\t' -1 1 -2 6 -o '1.1 2.1' <(cut -f1 mousePseudogenes.txt | sort -u) <(sort -t$'\t' -k6,6 $data/Ensembl/OrthoTables74/humanMouse.EnsemblOrtho.biomart80.map) | sort -u) > Mouse.pseudogenes.proteinCodingTx.humanOrtho
cut -f1 Mouse.pseudogenes.proteinCodingTx.humanOrtho | sort | uniq -c | awk '{if($1>1) print $2}' | sort -u > nonUniqMouse
cut -f2 Mouse.pseudogenes.proteinCodingTx.humanOrtho | sort | uniq -c | awk '{if($1>1) print $2}' | sort -u > nonUniqHuman
grep -v -f nonUniqMouse Mouse.pseudogenes.proteinCodingTx.humanOrtho | grep -v -f nonUniqHuman > tmp; mv tmp Mouse.pseudogenes.proteinCodingTx.humanOrtho
rm -rf nonUniqMouse nonUniqHuman

# Step 2 - Expand human genes with transcripts and gene symbols, then remove incomplete transcripts
join -t$'\t' -1 2 -2 2 -o '1.1 1.2 1.3 2.1' <(join -t$'\t' -1 2 -2 1 -o '1.1 2.1 2.2' <(sort -t$'\t' -k2 Mouse.pseudogenes.proteinCodingTx.humanOrtho) <(sort $data/Ensembl/geneTranscript.ensembl74.tsv) | sort -u | sort -t$'\t' -k2,2) <(sort -t$'\t' -k2,2 $data/Ensembl/humanGeneSymbol.humanEnsembl.biomart74.NoSyn.map) | sort -u > tmp
join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4' <(hgsql hg19 -Ne "SELECT name FROM ensGene WHERE cdsStartStat='cmpl' AND cdsEndStat='cmpl'" | sort) <(sort -t$'\t' -k3,3 tmp) > Mouse.pseudogenes.proteinCodingTx.humanOrtho 
rm -rf tmp

# Step 3 - Expand with the features
# chains (2.3), stop-codons (2.4), frameshifts (2.11), number of deleted exons (2.5), number of total exons (2.15), number of deleted bases (2.10), number of deleted bases (2.14), percent.id (2.16)
# then, expand with the KaKs values
join -t$'\t' -1 3 -2 2 -o '2.1 1.1 1.2 1.3 1.4 2.3 2.4 2.11 2.5 2.15 2.10 2.14 2.16' <(sort -t$'\t' -k3,3 Mouse.pseudogenes.proteinCodingTx.humanOrtho) <(grep ${assembly} $data/hg19.calltable.unique.placentals | sort -t$'\t' -k2,2) | sort -t$'\t' -k4,4 > tmp
join -t$'\t' -1 4 -2 2 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 2.2 2.3 2.11' <(sort -t$'\t' -k4,4 tmp) <(sort -t$'\t' -k2,2 $data/KaKs/speciesData/transcriptNonSynRate.${assembly}.tsv) | awk -F'\t' '{if($6==$15) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9"\t"$10"\t"$9/$10"\t"$11"\t"$12"\t"$11/$12"\t"$13"\t"$16}' | sort -u | sort -k2,2 -k3,3 -k4,4 > Mouse.pseudogenes.proteinCodingTx.humanOrtho
rm -rf tmp

# Step 4 - run the analysis of the intact genes:
sed -e 's/Rscript/#Rscript/g' $GENELOSS/src/sh/intactGeneTables.sh > intactGeneTables.sh
chmod 755 intactGeneTables.sh
./intactGeneTables.sh Mouse

rm -rf tmp
