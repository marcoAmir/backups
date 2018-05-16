#!/bin/bash -e

# this script check that fastq1 and fastq2 from paired-end library are synced. i.e., mate IDs are synced and sorted
# IMPORTANT: Note that if the fastq files need some processing (removal of disconcordant reads), the files will be overrided

fastq1=$1
fastq2=$2
base=`echo $fastq1 | awk -F'/' '{print $NF}' | cut -d"." -f1 | cut -d"_" -f1`

zcat ${fastq1} | grep "@${base}" | cut -d" " -f3 | sort | uniq -c | awk -v b=${base} '{print "\t"b"_1.fastq.gz: "$1"\treads;\t"$2}' > tmp1
zcat ${fastq2} | grep "@${base}" | cut -d" " -f3 | sort | uniq -c | awk -v b=${base} '{print "\t"b"_2.fastq.gz: "$1"\treads;\t"$2}' > tmp2

cat tmp1 tmp2

n1=`awk '{print $2}' tmp1`
n2=`awk '{print $2}' tmp2`

if [ $n1 -ne $n2 ];then
	echo -e "\t...number of reads is disconcordant: fixing FastQ mate pairs...\n"
	zcat ${fastq1} | awk '{printf ("%s",$0); n++; if(n%4==0){print ""} else {printf "\t"}}' | sort -u > sorted1.txt
	zcat ${fastq2} | awk '{printf ("%s",$0); n++; if(n%4==0){print ""} else {printf "\t"}}' | sort -u > sorted2.txt
	join -t$'\t' -1 1 -2 1 -o '1.1 1.2 1.3 1.4 2.1 2.2 2.3 2.4' sorted1.txt sorted2.txt > joined_tmp
	awk -F'\t' '{print $1; print $2; print $3; print $4}' joined_tmp > tmp1_fastq
	awk -F'\t' '{print $5; print $6; print $7; print $8}' joined_tmp > tmp2_fastq
	gzip tmp1_fastq tmp2_fastq
	mv tmp1_fastq.gz $fastq1
	mv tmp2_fastq.gz $fastq2
	zcat ${fastq1} | grep "@${base}" | cut -d" " -f3 | sort | uniq -c | awk -v b=${base} '{print "\t"b"_1.fastq.gz: "$1"\treads;\t"$2}' > tmp1
	zcat ${fastq2} | grep "@${base}" | cut -d" " -f3 | sort | uniq -c | awk -v b=${base} '{print "\t"b"_2.fastq.gz: "$1"\treads;\t"$2}' > tmp2
	cat tmp1 tmp2
	rm -rf sorted1.txt sorted2.txt joined_tmp
fi

rm -rf tmp1 tmp2 

m=`diff <(zcat ${fastq1} | grep "@${base}" | cut -d" " -f1) <(zcat ${fastq2} | grep "@${base}" | cut -d" " -f1) | wc -l`

if [ $m -eq 0 ];then
	echo -e "\t...mate IDs in ${fastq1} and ${fastq2} are consistent...\n"
else
	m=`diff <(zcat ${fastq1} | grep "@${base}" | cut -d" " -f1 | sort -u) <(zcat ${fastq2} | grep "@${base}" | cut -d" " -f1 | sort -u) | wc -l`
	if [ $m -eq 0 ];then
		echo -e "\t...Warning! mate IDs in ${fastq1} and ${fastq2} are NOT aligned...\n"
	else
		exit "\t...ERROR: mate IDs in ${fastq1} and ${fastq2} are NOT consistent...\n"
	fi
fi


