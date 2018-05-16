#!/bin/bash -e

if [ "$#" -ne 4 ] && [ "$#" -ne 3 ] && [ "$#" -ne 2 ]; then
	echo -e "\ncalling lastz one by one\n  Usage:\n  $0 target_assembly fasta_files_list\n"
	exit 1
fi

target=$1
in_fa_list=$2

# parameters
step=10
scores_file=loxAfr3_triMan1.scores


while read curr_fa
do
	base=`echo ${curr_fa} | sed -e "s/.fa//g"`
	echo -e "\n\t...processing ${base}"
	time lastz_32 /cluster/gbdb-bej/${target}/${target}.2bit[multiple,unmask] ${curr_fa} \
		--notransition --noytrim --step=${step} --scores=${scores_file} --format=sam \
		--output=${base}.sam
done < ${in_fa_list}
