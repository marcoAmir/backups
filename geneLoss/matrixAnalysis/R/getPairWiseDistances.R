require(ape)

hg19multiZ            <- read.tree('hg19.100way.nh')                                 # the 100-way multiz alignement newick
pairwiseBranchLengths <- cophenetic.phylo(hg19multiZ)

species   <- colnames(pairwiseBranchLengths)
human_ind <- which(species=="hg19")

for (i in 1:length(species)){
	x <- paste('hg19',species[i],pairwiseBranchLengths[human_ind, i], sep="\t")
	write(x, file='pairWiseDistances.tsv', append=TRUE)
}
