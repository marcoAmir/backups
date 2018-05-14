require(ggplot2)

source('/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/src/R/getClades.R')

geneFrameshift <- read.table('frameshifts.txt', col.names=c('ensgID','assembly','call'))
genes    <- as.character(unique(geneFrameshift$ensgID[which(geneFrameshift$call==1)]))	# extract the list of genes with detected losses

# For each gene I'll populate tree with calls for Lof (1) or otherwise (0)
# I'll iterate over each internal node and ask if the downstream clade contains only 1's and no zero - and if does I'll list the gene, and number of singel clade species that have lost it

queryTree$tip.label <- gsub("\\d","",queryTree$tip.label)	# strip all the assembly numbers from the species assembly abbreviations

genesOfMutualLoss <- c()	# this empty vector will contain the ensembl ID's of 
genesOfMultipleLoss <- c()	# this vector will contain all genes with multpile gene-losses
cladeMutation <- data.frame(ensgID=character(), clade=character())

c <- 1
for (currGene in genes){
	print(paste(c, currGene,sep='    '))
 	c <- c+1
	geneTree <- queryTree
	u <- rep(0, length(queryTree$tip.label))
	for (species in as.character(subset(geneFrameshift, ensgID==currGene & call==1)$assembly)){	# the list of species lost the particular gene	
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
			cladeMutation   <- rbind(cladeMutation, data.frame(ensgID=currGene, clade=clade))
		}
	}
}

# Now organize so you will have the most inclusive list of neighboring species sharing the same mutation
mutationSharedByClade <- data.frame(ensgID=character(), clade=character())
for (i in 1:nrow(cladeMutation)){ 
	print(i)
	mutation = as.character(cladeMutation$ensgID[i])
	clade = as.character(cladeMutation$clade[i])
	QA = 0	#QA
	if (mutation=='ENSG00000103310_chr16_21217087_21217088'){ #QA
		QA =1	#QA
		print(paste(mutation, clade, sep="    ")) #QA
	}	#QA
	if (mutation %in% mutationSharedByClade$ensgID){
		j <- which(mutationSharedByClade$ensgID==mutation)
		new_clade <- 1
		for (j_curr in j){
			u0 <- unique(unlist(strsplit(as.character(mutationSharedByClade$clade[j_curr]),'; ')))
			u1 <- unique(unlist(strsplit(as.character(cladeMutation$clade[i]),'; ')))
			if(length(u1) > length(u0) && length(intersect(u0, u1))){
				mutationSharedByClade$clade[j_curr] <- as.character(paste(u1, collapse='; '))
			}
			if(length(u1) == length(u0) && length(intersect(u0, u1) == length(u1))){
				next
			}
			if((!length(intersect(u0, u1)))){
				new_clade <- new_clade*1
			} else {
				new_clade <- new_clade*0
			}
		}
		if(new_clade && !(clade %in% unique(as.character(subset(mutationSharedByClade, ensgID==mutation)$clade)))){
			mutationSharedByClade <- rbind(mutationSharedByClade, data.frame(ensgID=mutation, clade=as.character(paste(u1, collapse='; '))))
		}
	} else {
		mutationSharedByClade <- rbind(mutationSharedByClade, data.frame(ensgID=mutation, clade=as.character(cladeMutation$clade[i])))
	}
	mutationSharedByClade <- unique(mutationSharedByClade)
	if (QA){	#QA
		print(subset(mutationSharedByClade, ensgID=='ENSG00000103310_chr16_21217087_21217088'))	#QA
	} #QA
}

write.table(mutationSharedByClade, file='cladeMutation.txt', row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# write doen genes that have multiple losses in the placental phylogeny
#for (g in setdiff(unique(genesOfMultipleLoss), unique(genesOfMutualLoss))){
#	write(g, file='multipleLossesGenes.txt', append=TRUE)
#}

# And, finally the script rewrites the input geneFrameshift file with the genes of mutual losses (and including the conserving species)


