require(ggplot2)
#library("scatterplot3d")
source('/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/src/R/multiplot.R')

species <- QUERY
geneTable <- read.table(paste0(species,'.intact.proteinCodingTx.humanOrtho'), col.names=c('assembly','qGeneID','humanGene','transcript','name','stop','fs','delExons','totalExons','fracDelExons','delCodingBases','totalCodingBases','fracDelBases','percentId','KaKs'))
geneTableRemain <- read.table(paste0(species,'.remainingCalls.proteinCodingTx.humanOrtho'), col.names=c('assembly','humanGene','transcript','name','stop','fs','delExons','totalExons','fracDelExons','delCodingBases','totalCodingBases','fracDelBases','percentId','KaKs'))
#geneTableLoss <- read.table(paste0(species,'.lossCalls.proteinCodingTx.humanOrtho'), col.names=c('assembly','humanGene','transcript','name','stop','fs','delExons','totalExons','fracDelExons','delCodingBases','totalCodingBases','fracDelBases','percentId','KaKs'))

curated <- read.table(paste0(species,'.curated'), col.names=c('curatedGenes'))

# geneTableCleanest - will hold a single transcript from each gene with the minimal perturbations
geneTableCleanest <- subset(geneTable, qGeneID %in% as.character(subset(data.frame(table(geneTable$qGeneID)), Freq==1)$Var1))	# First, I'll take all the one-to-one-to-one mappings
geneTableRemainCleanest <- subset(geneTableRemain, humanGene %in% as.character(subset(data.frame(table(geneTableRemain$humanGene)), Freq==1)$Var1))
#geneTableLossCleanest <- subset(geneTableLoss, humanGene %in% as.character(subset(data.frame(table(geneTableLoss$humanGene)), Freq==1)$Var1))

qGenes <- setdiff(as.character(unique(geneTable$qGeneID)), as.character(unique(geneTableCleanest$qGeneID)))
for (currGene in qGenes){
	print(currGene)
	currMap <- subset(geneTable, qGeneID==currGene)
	currMap <- currMap[order(-currMap$totalCodingBases),]
	geneTableCleanest <- rbind(geneTableCleanest, currMap[1,])
	#if(nrow(currMap)>1){
	#	geneTableCleanest <- rbind(geneTableCleanest, currMap[2,])
	#	next
	#}
}

hGenes <- setdiff(as.character(unique(geneTableRemain$humanGene)), as.character(unique(geneTableRemainCleanest$humanGene)))
for (currGene in hGenes){
	print(currGene)
	currMap <- subset(geneTableRemain, humanGene==currGene)
	currMap <- currMap[order(-currMap$totalCodingBases),]
	geneTableRemainCleanest <- rbind(geneTableRemainCleanest, currMap[1,])
}

#lGenes <- setdiff(as.character(unique(geneTableLoss$humanGene)), as.character(unique(geneTableLossCleanest$humanGene)))
#for (currGene in lGenes){
#	print(currGene)
#	currMap <- subset(geneTableLoss, humanGene==currGene)
#	currMap <- currMap[order(-currMap$totalCodingBases),]
#	geneTableLossCleanest <- rbind(geneTableLossCleanest, currMap[1,])
#}

geneTableCurated <- rbind(subset(geneTableRemain, humanGene %in% as.character(unique(curated$curatedGenes))), subset(geneTable[,-2], humanGene %in% as.character(unique(curated$curatedGenes))))

# Plotting the intact genes
p1 <- ggplot() + geom_point(data=geneTable, aes(x=jitter(stop), y=jitter(fs)), alpha=0.4) + geom_point(data=geneTableCleanest, aes(x=jitter(stop), y=jitter(fs)), color='steelblue', alpha=0.2) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Premature stop codons') + ylab('Frameshifts')
p2 <- ggplot() + geom_point(data=geneTable, aes(x=jitter(stop), y=fracDelBases), alpha=0.4) + geom_point(data=geneTableCleanest, aes(x=jitter(stop), y=fracDelBases), color='steelblue', alpha=0.2) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Premature stop codons') + ylab('fraction of deleted bases')
p3 <- ggplot() + geom_point(data=geneTable, aes(x=jitter(fs), y=fracDelBases), alpha=0.4) + geom_point(data=geneTableCleanest, aes(x=jitter(fs), y=fracDelBases), color='steelblue', alpha=0.2) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Frameshifts') + ylab('fraction of deleted bases')
p4 <- ggplot() + geom_point(data=geneTable, aes(x=jitter(fs), y=KaKs), alpha=0.4) + geom_point(data=geneTableCleanest, aes(x=jitter(fs), y=KaKs), color='steelblue', alpha=0.2) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Frameshifts') + ylab('Ka/Ks')
p5 <- ggplot() + geom_point(data=geneTable, aes(x=jitter(stop), y=KaKs), alpha=0.4) + geom_point(data=geneTableCleanest, aes(x=jitter(stop), y=KaKs), color='steelblue', alpha=0.2) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Premature stop codons') + ylab('Ka/Ks')
p6 <- ggplot() + geom_point(data=geneTable, aes(x=fracDelBases, y=KaKs), alpha=0.4) + geom_point(data=geneTableCleanest, aes(x=fracDelBases, y=KaKs), color='steelblue', alpha=0.2) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('fraction of deleted bases') + ylab('Ka/Ks')
png(paste0('intactGene.',species,'.png'))
multiplot(p1, p3, p2, p4, p5, p6, cols=3)
dev.off()

# Plotting the intact genes (longest transcript) on top of the transcripts that lack orthology to intact genes
p1 <- ggplot() + geom_point(data=geneTableRemain, aes(x=jitter(stop), y=jitter(fs)), color='red', alpha=0.7) + geom_point(data=geneTable, aes(x=jitter(stop), y=jitter(fs)), color='steelblue', alpha=0.3) + geom_point(data=geneTableCurated, aes(x=jitter(stop,0.1), y=jitter(fs,0.1)), color='green3', size=3, alpha=0.7) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Premature stop codons') + ylab('Frameshifts')
p2 <- ggplot() + geom_point(data=geneTableRemain, aes(x=jitter(stop), y=fracDelBases), color='red', alpha=0.7) + geom_point(data=geneTable, aes(x=jitter(stop), y=fracDelBases), color='steelblue', alpha=0.3) + geom_point(data=geneTableCurated, aes(x=jitter(stop,0.1), y=fracDelBases), color='green3', size=3, alpha=0.7) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Premature stop codons') + ylab('fraction of deleted bases')
p3 <- ggplot() + geom_point(data=geneTableRemain, aes(x=jitter(fs), y=fracDelBases), color='red', alpha=0.7) + geom_point(data=geneTable, aes(x=jitter(fs), y=fracDelBases), color='steelblue', alpha=0.3) + geom_point(data=geneTableCurated, aes(x=jitter(fs,0.1), y=fracDelBases), color='green3', size=3, alpha=0.7) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Frameshifts') + ylab('fraction of deleted bases')
p4 <- ggplot() + geom_point(data=geneTableRemain, aes(x=jitter(fs), y=KaKs), color='red', alpha=0.7) + geom_point(data=geneTable, aes(x=jitter(fs), y=KaKs), color='steelblue', alpha=0.3) + geom_point(data=geneTableCurated, aes(x=jitter(fs,0.1), y=KaKs), color='green3', size=3, alpha=0.7) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Frameshifts') + ylab('Ka/Ks')
p5 <- ggplot() + geom_point(data=geneTableRemain, aes(x=jitter(stop), y=KaKs), color='red', alpha=0.7) + geom_point(data=geneTable, aes(x=jitter(stop), y=KaKs), color='steelblue', alpha=0.3) + geom_point(data=geneTableCurated, aes(x=jitter(stop,0.1), y=KaKs), color='green3', size=3, alpha=0.7) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Premature stop codons') + ylab('Ka/Ks')
p6 <- ggplot() + geom_point(data=geneTableRemain, aes(x=fracDelBases, y=KaKs), color='red', alpha=0.7) + geom_point(data=geneTable, aes(x=fracDelBases, y=KaKs), color='steelblue', alpha=0.3) + geom_point(data=geneTableCurated, aes(x=fracDelBases, y=KaKs), color='green3', size=3, alpha=0.7) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('fraction of deleted bases') + ylab('Ka/Ks')
png(paste0('intactAndNonOrthoGenes.',species,'.png'))
multiplot(p1, p3, p2, p4, p5, p6, cols=3)
dev.off()

p1log <- p1 + theme(axis.text.x = element_text(angle=90, vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + scale_y_log10() + scale_x_log10()                    
p2log <- p2 + theme(axis.text.x = element_text(angle=90, vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + scale_x_log10()
p3log <- p3 + theme(axis.text.x = element_text(angle=90, vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + scale_x_log10()
p4log <- p4 + theme(axis.text.x = element_text(angle=90, vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + scale_x_log10()
p5log <- p5 + theme(axis.text.x = element_text(angle=90, vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + scale_x_log10()
p6log <- p6 + theme(axis.text.x = element_text(angle=90, vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) 
png(paste0('intactAndNonOrthoGenes.log.',species,'.png'))
multiplot(p1log, p3log, p2log, p4log, p5log, p6log, cols=3)
dev.off()

pdf(paste0('intactAndNonOrthoGenes.',species,'.pdf'), onefile = TRUE)
p1
p1log
p2
p2log
p3
p3log
p4
p4log
p5
p5log
p6
dev.off()

# Plotting the intact genes (longest transcript) on top of the transcripts that lack orthology to intact genes
#p1 <- ggplot() + geom_point(data=geneTableRemainCleanest, aes(x=jitter(stop), y=jitter(fs)), color='red', alpha=0.7) + geom_point(data=geneTableCleanest, aes(x=jitter(stop), y=jitter(fs)), color='steelblue', alpha=0.3) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Premature stop codons') + ylab('Frameshifts')
#p2 <- ggplot() + geom_point(data=geneTableRemainCleanest, aes(x=jitter(stop), y=fracDelBases), color='red', alpha=0.7) + geom_point(data=geneTableCleanest, aes(x=jitter(stop), y=fracDelBases), color='steelblue', alpha=0.3) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Premature stop codons') + ylab('fraction of deleted bases')
#p3 <- ggplot() + geom_point(data=geneTableRemainCleanest, aes(x=jitter(fs), y=fracDelBases), color='red', alpha=0.7) + geom_point(data=geneTableCleanest, aes(x=jitter(fs), y=fracDelBases), color='steelblue', alpha=0.3) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Frameshifts') + ylab('fraction of deleted bases')
#p4 <- ggplot() + geom_point(data=geneTableRemainCleanest, aes(x=jitter(fs), y=KaKs), color='red', alpha=0.7) + geom_point(data=geneTableCleanest, aes(x=jitter(fs), y=KaKs), color='steelblue', alpha=0.3) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Frameshifts') + ylab('Ka/Ks')
#p5 <- ggplot() + geom_point(data=geneTableRemainCleanest, aes(x=jitter(stop), y=KaKs), color='red', alpha=0.7) + geom_point(data=geneTableCleanest, aes(x=jitter(stop), y=KaKs), color='steelblue', alpha=0.3) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Premature stop codons') + ylab('Ka/Ks')
#p6 <- ggplot() + geom_point(data=geneTableRemainCleanest, aes(x=fracDelBases, y=KaKs), color='red', alpha=0.7) + geom_point(data=geneTableCleanest, aes(x=fracDelBases, y=KaKs), color='steelblue', alpha=0.3) + theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(angle=90, vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('fraction of deleted bases') + ylab('Ka/Ks')
#png(paste0('intactAndNonOrthoGenes.',species,'.png'))
#multiplot(p1, p3, p2, p4, p5, p6, cols=3)
#dev.off()

#pcaMatrix <- cbind(geneTableCleanest[,c(3,4,5,6,7,10,11)], type=rep('intact',nrow(geneTableCleanest)))
#pcaMatrix <- rbind(pcaMatrix, cbind(geneTableRemainCleanest[,c(1,2,3,4,5,8,9)], type=rep('nonOrtho', nrow(geneTableRemainCleanest))))
#pcaMatrix <- pcaMatrix[which(!is.na(pcaMatrix$KaKs)),]


#for (currGene in qGenes){
#	print(currGene)
#	currMap <- subset(geneTable, qGeneID==currGene)
#	currMap1st <- subset(currMap, stop==0 & fs==0 & fracDelBases==0)
#	if (nrow(currMap1st)==1){
#		geneTableCleanest <- rbind(geneTableCleanest, currMap1st)
#		c <- c+1
#		next
#	}
#	if (nrow(currMap1st)>1){
#		currMap1st <- currMap1st[order(-currMap1st$totalCodingBases, currMap1st$KaKs),]
#		geneTableCleanest <- rbind(geneTableCleanest, currMap1st[1,])
#		c <- c+1
#		next
#	}
#}
