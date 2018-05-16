#!/usr/bin/python

# Sept 2017 - this script extracts sequencing reads by their ID from an NGS project

import sys
import gzip

input_fq = sys.argv[1]
id_file  = sys.argv[2]	# file of read IDs to extract

input_name = input_fq.split('/')[-1].replace('.fastq.gz','.')

#mate = ''
#if '_1.fastq' in input_fq:
#	mate = '_1'
#if '_2.fastq' in input_fq:
#	mate = '_2'
#output_fq_file = input_name + id_file.replace('.baited','') + mate + '.fastq'
#output_fa_file = input_name + id_file.replace('.baited','') + mate + '.fa'

output_fq_file =  ('/').join(id_file.split('/')[:-1])
output_fa_file =  ('/').join(id_file.split('/')[:-1]) 
if len(id_file.split('/'))>1:
	output_fq_file =  output_fq_file + '/' + input_name + id_file.split('/')[-1].replace('.baited','') + '.fastq'
	output_fa_file =  output_fa_file + '/' + input_name + id_file.split('/')[-1].replace('.baited','') + '.fa'
else:
	output_fq_file =  input_name + id_file.split('/')[-1].replace('.baited','') + '.fastq'
	output_fa_file =  input_name + id_file.split('/')[-1].replace('.baited','') + '.fa'
	
wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))

out_fq = open(output_fq_file, 'w')
out_fa = open(output_fa_file, 'w')

with gzip.open(input_fq) as fq:
	line = fq.readline()
	l = 0
	while line:
		print l
		if l%4==0:
			data = line.strip().split()
			if data[0][1:] in wanted:
				out_fq.write(line)
				out_fa.write(data[0].replace('@', '> ') + '\n')
				line = fq.readline()
				l += 1
				out_fq.write(line)
				out_fa.write(line)
				for i in range(2):
					line = fq.readline()
					l += 1
					out_fq.write(line)
		line = fq.readline()
		l += 1
		
out_fq.close()
out_fa.close()
