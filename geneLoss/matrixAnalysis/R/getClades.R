require(ggplot2)
require(ape)

# This script will list all monophyletic groups in the placental mammals substree from the 100-way multiZ alignment

placentalLast <- 'dasNov3'

#queryTree <- read.tree('/cluster/u/amirma/geneLoss/data/UCSC/multiz100way/hg19.100way.nh')
queryTree <- read.tree('/cluster/u/amirma/geneLoss/data/UCSC/multiz100way/hg19.100way.topology.nh')

# Extract the placental mammal subtree
for(i in (length(queryTree$tip.label)+1):max(queryTree$edge)){
	queryTreeSub <- extract.clade(queryTree, i)
	if((placentalLast %in% queryTreeSub$tip.label) & ('hg19' %in% queryTreeSub$tip.label) & !('monDom5' %in% queryTreeSub$tip.label)){
		break
	}
}
queryTree <- queryTreeSub

nLeafs    <- length(queryTree$tip.label)

# Now iterate over each internal node to find subclades with 2, 3, 4, 5... species and to list them
for (cladeOrder in 2:40){
	print(paste("Clade order: ",cladeOrder,sep=""))
	for(i in (length(queryTree$tip.label)+1):max(queryTree$edge)){
		queryTreeSub <- extract.clade(queryTree, i)
		if(length(queryTreeSub$tip.label)==cladeOrder){
			write(c(queryTreeSub$tip.label), file='clades.txt', ncolumns=cladeOrder, append=TRUE, sep="\t")
		}
	}
}
