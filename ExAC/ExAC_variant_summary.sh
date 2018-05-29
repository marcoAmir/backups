#!/bin/bash -e 

# This script gets a genomic region (chr, start, end) and ExAC variants (vcf).
#  It will parse the info column to give a summary
#  that is identical to the summary given the ExAC browser (population frequency, etc.)

chr=$1
startPos=$2
endPos=$2
inVCF=$4

tabix ${inVCF} ${chr}:${startPos}-${endPos} > tmp.vcf

#./ExAC_variant_summary.py tmp.vcf




