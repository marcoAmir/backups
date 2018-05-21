# Example of obtaining sra data from ncbi and aligning it to a reference genome with Bowtie2 

1) Project: Duplex sequencing of Cockayne syndrome type B fibroblasts treated w/0J UVC 
	https://www.ncbi.nlm.nih.gov/sra/SRX4089536[accn]
	Name: CSB(1) 0J/m2
	Instrument: Illumina HiSeq 2000
	Strategy: OTHER
	Source: GENOMIC
	Selection: Hybrid Selection
	Layout: PAIRED

2) Get accession list (SRR_Acc_List.txt) of runs at:
	https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP147558
	(just analyze some of the samples - comment out experiments samples you don't want)
   And download the data:
	

