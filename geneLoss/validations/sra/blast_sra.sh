#!/bin/bash -e 

# Use this script to blast against an SRA file(s) with a given bait (fasta format)

sra=$1
bait_file=$2		# make sure you give absolute path
workdir=${PWD}

out_base=`echo ${bait_file} | awk -F'/' '{print $NF}' | sed -e "s/.fa//g"`
rm -rf ${out_base}*

if [ "$#" -ne 2 ]; then
  echo -e "\nBlast an SRA file with blastn_vdb\n  Usage:\n\t$0 sra_id bait_file(fa)\n(sra_id could be an sra file name an SRA id or list of SRA accesion - in the latter make sure it has txt suffix)"
  exit 1
fi

is_acc_file=`echo $sra | grep txt | wc -l`

if [ ${is_acc_file} -ne 0 ]; then
	while read sra_file
	do
		echo ${sra_file}
		base=`echo $sra_file | cut -d"." -f1`
		curr_dir=`echo ${PWD} | awk -F'/' '{print $NF}'`
		if [ $curr_dir != $base ]; then
			cd ${base}
		fi
		blastn_vdb -db ${sra_file}.sra -query ${bait_file} -out tmp.baited
		grep "^SRA:" tmp.baited | cut -d" " -f1 | sed -e "s/SRA://g" | sort -u >> ${workdir}/${out_base}.${sra_file}.baited
		cd ${workdir}		
	done < ${sra}
else
	is_sra_file=`echo $sra | grep sra | wc -l`
	base=`echo $sra | cut -d"." -f1`
	if [ ${is_sra_file} -ne 0 ]; then
		blastn_vdb -db ${sra} -query ${bait_file} -out tmp.baited 
		grep "^SRA:" tmp.baited | cut -d" " -f1 | sed -e "s/SRA://g" | sort -u >> ${workdir}/${out_base}.${base}.baited
	else
		curr_dir=`echo ${PWD} | awk -F'/' '{print $NF}'`
		if [ $curr_dir != $base ]; then
			cd ${base}
		fi
		blastn_vdb -db ${base}.sra -query ${bait_file} -out tmp.baited 
		grep "^SRA:" tmp.baited | cut -d" " -f1 | sed -e "s/SRA://g" | sort -u >> ${workdir}/${out_base}.${base}.baited
		cd ${workdir}
	fi
fi
