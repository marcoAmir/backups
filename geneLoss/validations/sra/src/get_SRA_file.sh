#!/bin/bash -e

# from the SRA factsheet (available as pdf at www.ncbi.nlm.nih.gov/sra/docs):
#   Given an XRR (SRR/ERR/DRR) accession, you can use the following steps to reconstruct the FTP path for the .sra file:
#    The base FTP path is ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/
#    Append /XRR to get to the different source directory (with X being S, E, or D)
#    Append /XRR### with the # being the first three digits of the XRR accession, for SRR1427233, use /SRR142
#    Append XRR full accession, for SRR1427233, use /SRR1427233
#    Append the full accession with .sra extension, for SRR1427233, use /SRR1427233.sra to arrive at:
#   ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR142/SRR1427233/SRR1427233.sra
#   For ascp, replace the ftp.ncbi.nlm.nih.gov with anonftp@ftp-private.ncbi.nlm.nih.gov: to arrive at:
#   anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR142/SRR1427233/SRR1427233.sra 

query=$1	# SRA accession. e.g., SRR1427233

srr_list=`echo $query | grep ".txt" | wc -l`

if [ $srr_list -ne 0 ];then
	while read curr_sra
	do
		mkdir ${curr_sra}
		cd ${curr_sra}
		base=`echo ${curr_sra} | sed -e "s/[0-9]//g"`
		subdir=`echo ${curr_sra} | cut -c1-6`
		ftp_path='ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/'${base}'/'${subdir}'/'${curr_sra}'/'${curr_sra}'.sra'
		wget ${ftp_path}
		cd ../
	done < $query
else
	mkdir ${query}
	cd ${query}
	base=`echo $query | sed -e "s/[0-9]//g"`
	subdir=`echo $query | cut -c1-6`
	ftp_path='ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/'${base}'/'${subdir}'/'${query}'/'${query}'.sra'
	wget ${ftp_path}
	cd ../
fi


