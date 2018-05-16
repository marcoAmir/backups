#!/bin/bash -e

# this script analyzes sra file: it extracts fastq, aligns to a target assembly with bowtie2, and get a bed file
# currently works only on paired-end illumina libraries

BT2_HOME=/cluster/u/amirma/bin/Bowtie2/bowtie2-2.3.2

####### Settings: ###########################################################################
torch_sra=0
torch_fastq=0
torch_sam=1
# The above flags, if non-zero, will delete data after finish execution #
skip_FastQC=1
# Read handling
ReadLengthTrim=150	# some libraries comes with varying read length 
				# whose quality decreases as the read extends. 
				# set a number >0 to handle your fastqs and trim the reads
				# if <=0 no trimming
mapQuality=10		# exclude uncertain maps from sam (in the conversion to bam and bed)
local_align=1		# if selected run local alignment
#############################################################################################

workdir=${PWD}

sra_file=$1				# input name for sra
indexed_genome=$2			# input target genome, indexed (e.g., /cluster/u/amirma/geneLoss/hg38/validations/sra/cow/PRJEB14827/indexed_genomes/bosTau8)
library=$3				# gets "single" or "paired"

if [ "$#" -ne 3 ] && [ "$#" -ne 2 ]; then
  echo -e "\nAnalysis of a convergent evolution experiment and its controls\n  Usage:\n  $0 sra_file indexed_genome single/paired(OPT)\n\t\tdefault: paired"
  exit 1
fi

if [ "$#" -eq 2 ]; then
	library=paired
fi
if [ ${library} != "paired" ] && [ ${library} != "single" ]; then
	library=paired
fi
echo ${library}

base=`echo $sra_file | cut -d"." -f1`
target=`echo $indexed_genome | awk -F'/' '{print $NF}'`

curr_dir=`echo ${PWD} | awk -F'/' '{print $NF}'`

if [ $curr_dir != $base ]; then
	cd ${base}
fi

# get the fastq files:
if [ ! -e *fastq.gz ]; then
	echo -e "\n\n\n\t...getting fastq files...\n"
	if [ ${library} = "single" ]; then
		fastq-dump --gzip ${sra_file}
	else
		fastq-dump --split-files --gzip ${sra_file}
	fi
fi

# trim reads if asked to (ReadLengthTrim>0)
if [ ${ReadLengthTrim} -gt 0 ]; then
	echo -e "\n\n\n\t...trimming read length to ${ReadLengthTrim} (WARNING: This will permanently modify your Fastqs)...\n"
	for f in *.fastq.gz
	do
		zcat -f ${f} | sed -e "s/length=[0-9]*/length=${ReadLengthTrim}/g" | awk -v L=${ReadLengthTrim} '{if(NR%2==0) {print substr($0,1,L)} else {print $0}}' > tmp; gzip tmp; mv tmp.gz ${f}
	done
fi

# checking paired-end fastq concordance
if [ ${library} = "paired" ]; then
	echo -e "\n\n\n\t...checking paired-end fastq concordance...\n"
	/cluster/u/amirma/geneLoss/hg38/validations/sra/src/check_pairedEndFastQs.sh ${PWD}/${base}_1.fastq.gz ${PWD}/${base}_2.fastq.gz
fi

# FastQC - quality check for the fastq files:
if [ $skip_FastQC -eq 0 ]; then
	echo -e "\n\n\n\t...quality check of fastq files (FastQC)...\n"
	mkdir QC
	for f in *.fastq.gz
	do
		fastqc ${f} -o QC
	done
fi

# Align to a selected genome with bowtie2
echo -e "\n\n\n\t...aligning reads to ${target}...\n"
echo ${indexed_genome}*bt2 | awk '{for(i=1; i<=NF; i++) print "ln -f -s "$i}' | bash
if [ ${library} = "paired" ]; then
	if [ $local_align -eq 0 ];then
		${BT2_HOME}/bowtie2 -x ${target} -p 4 -1 ${base}_1.fastq.gz -2 ${base}_2.fastq.gz -S ${base}_${target}.sam
	else
		${BT2_HOME}/bowtie2 --local -x ${target} -p 4 -1 ${base}_1.fastq.gz -2 ${base}_2.fastq.gz -S ${base}_${target}.sam
	fi
else
	if [ $local_align -eq 0 ];then
		${BT2_HOME}/bowtie2 -x ${target} -p 4 -U ${base}.fastq.gz -S ${base}_${target}.sam
	else
		${BT2_HOME}/bowtie2 --local -x ${target} -p 4 -U ${base}.fastq.gz -S ${base}_${target}.sam
	fi
fi
#${BT2_HOME}/bowtie2 -x ${target} -p 2 -1 ${base}_1.fastq.gz -2 ${base}_2.fastq.gz -S ${base}_${target}.sam

# Convert sam to bam, bam to bed
echo -e "\n\n\n\t...converting sam to bam...\n"
samtools view -bS -q ${mapQuality} ${base}_${target}.sam > ${base}_${target}.bam
echo -e "\n\n\n\t...converting bam to bed(12)...\n"
bedtools bamtobed -bed12 -i ${base}_${target}.bam > ${base}_${target}.bed
cat <(echo "track name=${base}_${target}_mapped_reads description=\"${base} reads\" color=160,160,160") <(cat ${base}_${target}.bed) > tmp; mv tmp ${base}_${target}.bed

# Clean up
if [ $torch_sra -ne 0 ];then
	rm -rf ${sra_file}
fi

if [ $torch_fastq -ne 0 ];then
	rm -rf *.fastq.gz
fi

if [ $torch_sam -ne 0 ];then
	rm -rf ${base}_${target}.sam
else
	gzip ${base}_${target}.sam
fi

gzip *.bed

rm -rf *.bt2

if [ $curr_dir != $base ];then
	cd ${workdir}
fi
	
