require(ggplot2)

source('/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/src/R/getClades.R')

geneLoss <- read.table('geneSpecies.placentals.call', col.names=c('ensgID','assembly','call'))
genes    <- as.character(unique(geneLoss$ensgID[which(geneLoss$call==1)]))	# extract the list of genes with detected losses

# For each gene I'll populate tree with calls for Lof (1) or otherwise (0)
# I'll iterate over each internal node and ask if the downstream clade contains only 1's and no zero - and if does I'll list the gene, and number of singel clade species that have lost it

queryTree$tip.label <- gsub("\\d","",queryTree$tip.label)	# strip all the assembly numbers from the species assembly abbreviations

genesOfMutualLoss <- c()	# this empty vector will contain the ensembl ID's of 
genesOfMultipleLoss <- c()	# this vector will contain all genes with multpile gene-losses
geneCladeLosses <- data.frame(ensgID=character(), clade=character())

c <- 1
for (currGene in genes){
	print(paste(c, currGene,sep='    '))
 	c <- c+1
	geneTree <- queryTree
	u <- rep(0, length(queryTree$tip.label))
	for (species in as.character(subset(geneLoss, ensgID==currGene & call==1)$assembly)){	# the list of species lost the particular gene	
		u[which(queryTree$tip.label==species)] <- 1
	}
	if(sum(u) > 1){
		genesOfMultipleLoss <- c(genesOfMultipleLoss, currGene)
	}
	geneTree$tip.label <- u
	# at this stage the tree for a single gene has it's leafs populated with 0's and 1's, and I'll now find:
	#	i - genes with internal nodes (clades) that has all leafs of 1 (Lof) 
	for(i in (length(geneTree$tip.label)+1):max(geneTree$edge)){
		geneTreeSub  <- extract.clade(geneTree, i)
		queryTreeSub <- extract.clade(queryTree, i)
		v <- as.numeric(geneTreeSub$tip.label)
		clade <- paste(queryTreeSub$tip.label, collapse="; ")
		if (length(unique(v))==1 && unique(v)==1){
			genesOfMutualLoss <- c(genesOfMutualLoss, currGene)
			geneCladeLosses   <- rbind(geneCladeLosses, data.frame(ensgID=currGene, clade=clade))
		}
	}
}

write.table(geneCladeLosses, file='geneCladeLosses.txt', row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# write doen genes that have multiple losses in the placental phylogeny
for (g in setdiff(unique(genesOfMultipleLoss), unique(genesOfMutualLoss))){
	write(g, file='multipleLossesGenes.txt', append=TRUE)
}

# And, finally the script rewrites the input geneLoss file with the genes of mutual losses (and including the conserving species)


