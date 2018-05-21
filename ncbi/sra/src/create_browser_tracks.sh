#!/bin/bash -e 

assembly=$1
chr=$2
startPos=$3
endPos=$4

outFile=${assembly}_${chr}
wigFile=${outFile}

prjID=`echo ${PWD} | awk -F'/' '{print $NF}'`

if [ "$#" -ne 4 ] && [ "$#" -ne 3 ] && [ "$#" -ne 2 ]; then
	echo -e "\nCreating a browser tracks from aligned reads\n  Usage:\n  $0 assembly targetChrom optional:startPos optional:endPos\n\t\tdefault: if startPos and/or endPos not selected will do the whole chromosome"
	exit 1
fi

if [ "$#" -eq 2 ] || [ "$#" -eq 3 ]; then
	startPos=0
	endPos=10000000000
	outFile=${outFile}.bed
	wigFile=${wigFile}.wig	
else
	outFile=${outFile}_${startPos}_${endPos}.bed
	wigFile=${wigFile}_${startPos}_${endPos}.wig
fi

rm -rf ${outFile}

while read currID
do
	echo -e "...${currID}..."
	currBed=${currID}/*${assembly}.bed*
	if [ $startPos -eq 0 ]; then
		zcat -f $currBed | grep ${chr} >> ${outFile}
	else	
		zcat -f $currBed | grep ${chr} | awk -v startPos=${startPos} -v endPos=${endPos} -F'\t' '{if($2>=startPos && $2<endPos) print $0}'>> ${outFile}
		zcat -f $currBed | grep ${chr} | awk -v startPos=${startPos} -v endPos=${endPos} -F'\t' '{if($3>startPos && $3<=endPos) print $0}'>> ${outFile}
		zcat -f $currBed | grep ${chr} | awk -v startPos=${startPos} -v endPos=${endPos} -F'\t' '{if($2<=startPos && $3>=endPos) print $0}'>> ${outFile}
		sort -u ${outFile} > tmp; mv tmp ${outFile}
	fi
done < SRR_Acc_List.txt

bedSort ${outFile} tmp.bed
/cluster/u/amirma/geneLoss/hg38/validations/sra/src/bed2wig.py -i tmp.bed --type bed -o tmp.wig
cat <(echo -e "track name=${prjID}_mapped_reads_${chr} description=\"${prjID} reads\" color=160,160,160") <(cat tmp.bed) > ${outFile}
cat <(echo -e "track name=${prjID}_mapped_reads_${chr}_pile description=\"${prjID} reads pile\" color=120,120,120") <(cat tmp.wig) > ${wigFile}
rm -rf tmp*
