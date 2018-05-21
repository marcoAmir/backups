#!/bin/bash -e

src=/cluster/u/amirma/geneLoss/hg38/validations/sra/src

project_list=$1		# a list of SRA runs (e.g., SRR_Acc_List.txt)
in_ids=$2		# list of read ids to extract

base=`echo $in_ids | sed -e "s/.baited//g"`

runs_dir=`echo $project_list | awk -F'/' '{for(i=1; i<NF; i++) printf $i"/"} END {print ""}'`

while read curr_sra
do
	echo -e "\n\t...extracting from ${runs_dir}/${curr_sra}/${curr_sra}_1.fastq.gz "
	${src}/extract_fqReads_by_ID.py  ${runs_dir}${curr_sra}/${curr_sra}_1.fastq.gz ${in_ids}
	echo -e "\n\t...extracting from ${runs_dir}/${curr_sra}/${curr_sra}_2.fastq.gz "
	${src}/extract_fqReads_by_ID.py  ${runs_dir}${curr_sra}/${curr_sra}_2.fastq.gz ${in_ids}
	cat *.fastq >> tmp.fastq;
	cat *.fa >> tmp.fa; 
done < ${project_list}

mv tmp.fastq ${base}.fastq
mv tmp.fa ${base}.fa
