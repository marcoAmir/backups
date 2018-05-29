#!/usr/bin/python

# Following documentation at:
#	http://pyvcf.readthedocs.io/en/latest/INTRO.html

import vcf

vcf_reader = vcf.Reader(open('MMP21.vcf', 'r'))			# Reader - for training purposes, the MMP21 entries from ExAC
#Metadata regarding the VCF file itself can be investigated through the following attributes:
#        Reader.metadata
#        Reader.infos
#        Reader.filters
#        Reader.formats
#        Reader.samples


record = next(vcf_reader)					# Record - the attributes of a Record are the 8 fixed fields from the VCF spec: Record.CHROM, Record.POS, Record.ID, Record.REF, Record.ALT, Record.QUAL. Record.FILTER, Record.INFO``

vcf_writer = vcf.Writer(open('MMP21_output.vcf', 'w'), vcf_reader)	# Writer - class providing a way of writing a VCF file. Currently, you must specify a template Reader (vcf_reader) which provides the metadata:

# Some convenient methods to examine properties for each Record (just some examples):
print record.is_snp, record.is_indel, record.is_transition, record.is_deletion	# returns a boolean
print record.var_type, record.var_subtype					# e.g., 'snp' or 'ts'
print record.is_monomorphic

# write a vcf:
for record in vcf_reader:
	vcf_writer.write_record(record)






