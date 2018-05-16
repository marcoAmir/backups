#!/usr/bin/python

from Bio import SeqIO
input_file = "big_file.fastq"
id_file = "short_list.txt"
output_file = "short_list.fastq"
wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))

print "Found %i unique identifiers in %s" % (len(wanted), id_file)
records = (r for r in SeqIO.parse(input_file, "fastq") if r.id in wanted)
count = SeqIO.write(records, output_file, "fastq")
print "Saved %i records from %s to %s" % (count, input_file, output_file)
if count < len(wanted):
	print "Warning %i IDs not found in %s" % (len(wanted)-count, input_file)
