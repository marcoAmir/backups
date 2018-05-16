#!/bin/bash -e 

in_list=$1

while read sam_file 
do
	echo -e "\n\t...processing ${sam_file}"
	bam_file=`echo $sam_file | sed -e "s/sam/bam/g"`
	bed_file=`echo $sam_file | sed -e "s/sam/bed/g"`
	samtools view -bS ${sam_file}> ${bam_file}
	bedtools bamtobed -bed12 -i ${bam_file} > ${bed_file}
	cut -f4 ${bed_file} | sort | uniq -c | awk '{if($1==1) print $2}' > uniquely_mapped_reads.tmp
	join -t$'\t' -1 4 -2 1 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12' \
		<(sort -t$'\t' -k4,4 ${bed_file}) <(sort -u uniquely_mapped_reads.tmp) > tmp ; mv tmp ${bed_file}
done < ${in_list}
