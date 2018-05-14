#!/bin/bash -e

# This script will take as an input a callTable file (24 column files found in this directory) and will parse it into a *.vec file that gives ternary annotations for each transcript.
# in the vec file, the first column is the transcript id, the second is a comma-delimited list of species in which the transcript is intact, the third is unknowns, and the fourth is species with LoF transcripts

if [ "$#" -ne 2 ]; then
  echo -e "\nUsage:\n  $0 callTableFile output\n"
  exit 1
fi

callTable=$1
output=$2

rm -rf allPairs.txt

tail -n+2 $callTable | cut -f2 | sort -u > transcripts.txt
tail -n+2 $callTable | cut -f1 | sort -u > assemblies.txt
while read a
do
	echo $a
	while read t
	do
		echo -e "$t\t$a" >> allPairs.txt
	done < transcripts.txt
done < assemblies.txt

# I - All the species-transcript pairs annotated as intacts (coloumn 24 is '0') goes to 'intacts'. 
#     Then, if there are duplicates with the unknowns - they will be marked as intacts (that means we find at least one chain in which the transcript is intact)
awk -F'\t' '{if($24==0) {sub(/[0-9]+/,"",$1); print $2"\t"$1}}' $callTable | sort -u | awk -F'\t' '{if(a[1]==$1) {printf ","$2} else {a[1]=$1; print ""; printf $0",hg"}}' | awk -F'\t' '{if($0!="") print $0}' > intacts
cat <(cat intacts) <(comm -3 -2 <(sort transcripts.txt) <(cut -f1 intacts | sort -u) | awk -F'\t' '{print $1"\thg"}') > tmp; mv tmp intacts

# II - All the unannotated species-transcript pairs and all the pairs with question marks goes to 'unknowns' 
comm -1 -3 <(awk -F'\t' '{sub("[0-9]+","",$1); print $2"\t"$1}' $callTable | sort -u) <(awk -F'\t' '{sub("[0-9]+","",$2); print $1"\t"$2}' allPairs.txt | sort -u) > unknowns
cat <(cat unknowns) <(awk -F'\t' '{if($24=="?") {sub(/[0-9]+/,"",$1); print $2"\t"$1}}' $callTable) | sort -u > tmp; mv tmp unknowns 
comm -1 -3 <(sed -e 's/,/\t/g' intacts | awk -F'\t' '{for(i=2; i<=NF; i++) print $1"\t"$i}' | sort -u) <(sort -u unknowns) > tmp; mv tmp unknowns
awk -F'\t' '{if(a[1]==$1) {printf ","$2} else {a[1]=$1; print ""; printf $0}}' unknowns | awk -F'\t' '{if($0!="") print $0}' > tmp; mv tmp unknowns
cat <(cat unknowns) <(comm -3 -2 <(sort transcripts.txt) <(cut -f1 unknowns | sort -u) | awk -F'\t' '{print $1"\t"}') > tmp; mv tmp unknowns

# III - all the species-transcript pairs that are annotated for Loss of function (Lof, column 24 is '1') goes to Lofs
# duplicates with unknowns and intact are removed
comm -3 -2 <(awk -F'\t' '{if($24==1) {sub(/[0-9]+/,"",$1); print $2"\t"$1}}' $callTable | sort) <(sort unknowns) | sort -u > Lofs
comm -1 -3 <(cat intacts unknowns | sed -e 's/,/\t/g' | awk -F'\t' '{for(i=2; i<=NF; i++) print $1"\t"$i}' | sort -u) <(sort -u Lofs) | sort -u > tmp; mv tmp Lofs
awk -F'\t' '{if(a[1]==$1) {printf ","$2} else {a[1]=$1; print ""; printf $0}}' Lofs | awk -F'\t' '{if($0!="") print $0}' > tmp; mv tmp Lofs
cat <(cat Lofs) <(comm -3 -2 <(sort transcripts.txt) <(cut -f1 Lofs | sort -u) | awk -F'\t' '{print $1"\t"}') > tmp; mv tmp Lofs

join -t$'\t' -1 1 -2 1 -o '1.1 1.2 1.3 2.2' <(join -t$'\t' -1 1 -2 1 -o '1.1 1.2 2.2' <(sort intacts) <(sort unknowns) | sort) <(sort Lofs) > $output

rm -rf transcripts.txt Lofs unknowns intacts allPairs.txt

