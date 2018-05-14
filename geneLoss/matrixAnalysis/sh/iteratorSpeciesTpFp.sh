#!/bin/bash -e


# This script will iterate over the ~30 placental mammals who has Ensembl orthology tables with human
# I - it will first generate geneSpecies call file (1 gene is lost, 0 gene is intact in a species) based on preference on the human transcripts, see:	/cluster/u/amirma/geneLoss/matrixAnalysis/placentalMammalLosses/humanFilterTypes.txt
# II - then, given preferences on orthologous species (see orthoFilterTypes.txt) it will partition each species into TP and FP, update the callsVsEnsembl file and generate scatter, and will generate a summary count for the FP calls in the data

# Inputs (order matters): 
# -----------------------  
# 	speciesAnalyzed.txt
#	human-filters	(as h1,h2)
#	ortho-filters	(as o1,o2)
#	frac_del_bases	(as float in range [0-1])

src=$GENELOSS/src

if [ "$#" -ne 4 ]; then
  echo -e "\nUsage:\n  $0 speciesFile hFilters oFilters delFrac\n"
  exit 1
fi

# I
speciesFile=$1
hfiltersApplied=$2
ofiltersApplied=$3
delFrac=$4

applyDelFrac=`echo $delFrac'>'0 | bc -l`

echo $hfiltersApplied
echo $ofiltersApplied

jobId=$RANDOM
outDir=${hfiltersApplied}_${ofiltersApplied}_${delFrac}_${jobId}
mkdir $outDir

h1=0	# requiring human ensembl transcript with complete state-end
h2=0	# requiring human ensembl transcript with complete state-end, with UTRs
h3=0	# filtering out any human transcripts that overlap with seg-dup regions
h4=0	# taking only the subset of the data, where the gene-losses are shared between neighboring taxa (clade losses)
h5=0	# all duplicated species-transcript anotations (because of calls from multiple chains) are treated as unknowns
h6=0	# Allow a more permisive definition of Lof at the transcript level
h7=0	# Exclude genes that were 'protein-coding' in Ensembl74 but are no longer in Ensembl81
o1=0	# requiring that the orthologous transcripts are complete
o2=0	# requiring that the orthologous transcripts are complete, and have UTRs

if [[ $hfiltersApplied == *"h1"* ]]
then
	h1=1
fi

if [[ $hfiltersApplied == *"h2"* ]]
then
	h1=0
	h2=1
fi

if [[ $hfiltersApplied == *"h3"* ]]
then
	h3=1
fi

if [[ $hfiltersApplied == *"h4"* ]]
then
	h4=1
fi

if [[ $hfiltersApplied == *"h5"* ]]
then
	h5=1
fi

if [[ $hfiltersApplied == *"h6"* ]]
then
	h6=1
fi

if [[ $hfiltersApplied == *"h7"* ]]
then
	h7=1
fi

echo -e "\n  ...converting transcript to gene-species calls..."
$src/sh/humanTranscriptsToGenelosses.sh $h1 $h2 $h3 $h4 $h5 $h6 $h7

#if [[ $h4 -ne 0 ]]
#then
#	Rscript $src/R/carveCladeLosses.R
#	join -t$'\t' -1 1 -2 1 -o '1.1 1.2 1.3' <(sort geneSpecies.placentals.call) <(sort geneCladeLosses.txt) | sort -u > tmp; mv tmp geneSpecies.placentals.call
#fi

if [[ $ofiltersApplied == *"o1"* ]]
then
	o1=1
fi

if [[ $ofiltersApplied == *"o2"* ]]
then
	o1=0
	o2=1
fi

# II
sed -e "s/hfiltersApplied='h0'/hfiltersApplied='$hfiltersApplied'/g;s/ofiltersApplied='h0'/ofiltersApplied='$hfiltersApplied'/g;s/frac=0.25/frac=$delFrac/g" $src/sh/featureAnalysis.sh > featureAnalysis.sh
chmod 755 featureAnalysis.sh

rm -rf TPs.txt FPs.txt totalCalls.txt
while read s
do 
	echo -e "\n$s"
	if [[ $s == *"#"* ]]    
	then
		continue
	fi
	species=`echo $s | cut -d" " -f2`
	assembly=`echo $s | cut -d" " -f1 | sed -e 's/[0-9]//g'`
	./featureAnalysis.sh $s $o1 $o2 $applyDelFrac
	t=`tail -n+2 $PWD/$species/TP/tpGeneTranscripts.tsv | wc -l`
	f=`tail -n+2 $PWD/$species/FP/fpGeneTranscripts.tsv | wc -l`
	tf=`awk -F'\t' '{if($3==1) print $1"\t"$2}' geneSpecies.placentals.call | grep -w $assembly | wc -l`
	if [[ $t -ne 0 ]]
	then
		tail -n+2 $PWD/$species/TP/tpGeneTranscripts.tsv | cut -f1,3 | sort -u >> TPs.txt
	fi
	if [[ $f -ne 0 ]]
	then
		tail -n+2 $PWD/$species/FP/fpGeneTranscripts.tsv | cut -f1,3 | sort -u >> FPs.txt
	fi
	if [[ $tf -ne 0 ]]
	then
		awk -F'\t' '{if($3==1) print $1"\t"$2}' geneSpecies.placentals.call | grep -w $assembly >> totalCalls.txt
	fi
	mv $species $outDir
done < $speciesFile

cp hg19.vec $outDir
cp geneSpecies.placentals.call $outDir

if [ -f geneCladeLosses.txt ];
then
	cp geneCladeLosses.txt $outDir
fi 

# III (update log):
if [ ! -e transToGene.filtersCounts.orthoFilters.log ]
then 
	echo -e "hfilters\tofilters\ttotalTranscripts\ttotalGenes\tGeneLosses\tPlacentalLosses\tnGenes\tnSpecies\tensPlacentalLosses\tTPensPlacentalLosses\tFPensPlacentalLosses\tThresfractionDeletedBases\tFDR\tJob" > transToGene.filtersCounts.orthoFilters.log
fi
transToGene=`tail -1 transToGene.filtersCounts.log | cut -f2-7`
ensPlacentalLosses=`wc -l totalCalls.txt | cut -d" " -f1`
TPensPlacentalLosses=`wc -l TPs.txt | cut -d" " -f1`
FPensPlacentalLosses=`wc -l FPs.txt | cut -d" " -f1`
ensPlacentalLosses2=$((TPensPlacentalLosses+FPensPlacentalLosses))
echo -e "$hfiltersApplied\t$ofiltersApplied\t$transToGene\t$ensPlacentalLosses2\t$TPensPlacentalLosses\t$FPensPlacentalLosses\t$delFrac" | awk -F'\t' '{print $0"\t"$11/$9"\t"'$jobId'}' >> transToGene.filtersCounts.orthoFilters.log

