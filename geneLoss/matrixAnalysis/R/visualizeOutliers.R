require(ggplot2)

source('/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/src/R/multiplot.R')

jitterFactor <- 3

gamma <- 0.1
nu <- 0.1

species <- 'Mouse'
geneTable <- read.table(paste0('/cluster/u/amirma/geneLoss/intactGenes/',species,'.intact.proteinCodingTx.canonical.humanOrtho'), col.names=c('assembly','qGeneID','humanGene','transcript','name','stop','fs','delExons','totalExons','fracDelExons','delCodingBases','totalCodingBases','fracDelBases','percentId','KaKs'))
geneTableRemain <- read.table(paste0('/cluster/u/amirma/geneLoss/intactGenes/',species,'.remainingCalls.proteinCodingTx.canonical.humanOrtho'), col.names=c('assembly','humanGene','transcript','name','stop','fs','delExons','totalExons','fracDelExons','delCodingBases','totalCodingBases','fracDelBases','percentId','KaKs'))

curated <- read.table(paste0('/cluster/u/amirma/geneLoss/intactGenes/',species,'.curated'), col.names=c('curatedGenes'))

geneTableCleanest <- subset(geneTable, qGeneID %in% as.character(subset(data.frame(table(geneTable$qGeneID)), Freq==1)$Var1))	# First, I'll take all the one-to-one-to-one mappings
geneTableRemainCleanest <- subset(geneTableRemain, humanGene %in% as.character(subset(data.frame(table(geneTableRemain$humanGene)), Freq==1)$Var1))

qGenes <- setdiff(as.character(unique(geneTable$qGeneID)), as.character(unique(geneTableCleanest$qGeneID)))
for (currGene in qGenes){
	print(currGene)
	currMap <- subset(geneTable, qGeneID==currGene)
	currMap <- currMap[order(-currMap$totalCodingBases),]
	geneTableCleanest <- rbind(geneTableCleanest, currMap[1,])
}

hGenes <- setdiff(as.character(unique(geneTableRemain$humanGene)), as.character(unique(geneTableRemainCleanest$humanGene)))
for (currGene in hGenes){
	print(currGene)
	currMap <- subset(geneTableRemain, humanGene==currGene)
	currMap <- currMap[order(-currMap$totalCodingBases),]
	geneTableRemainCleanest <- rbind(geneTableRemainCleanest, currMap[1,])
}

geneTableCurated <- rbind(subset(geneTableRemain, humanGene %in% as.character(unique(curated$curatedGenes))), subset(geneTable[,-2], humanGene %in% as.character(unique(curated$curatedGenes))))


###############################################################
# Now read the results from the outlier prediction:
outliersPred <- read.table(paste0('output/TxPredictions.',species,'.One-Class_SVM.gamma_',gamma,'.nu_',nu), col.names=c('assembly','humanGene','transcript','name','outlier'))
outlierTXs   <- as.character(subset(outliersPred, outlier==1)$transcript)  
geneTableOutliers        <- geneTable[geneTable$transcript %in% outlierTXs,]
geneTableRemainOutliers  <- geneTableRemain[geneTableRemain$transcript %in% outlierTXs,]
geneTableCuratedOutliers <- geneTableCurated[geneTableCurated$transcript %in% outlierTXs,]
sampOutliers <- sample(which(geneTableOutliers$fs+geneTableOutliers$stop<=1), 10)	# sample outliers from the 'intact' genes (row indices in geneTableOutliers)
###############################################################

# frameshifts Vs. stop codons
p1 <- ggplot() + 
	geom_point(data=geneTableRemain, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor)), alpha=0.3) + 
	geom_point(data=geneTable, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor)), color='steelblue', alpha=0.3) + 
	geom_point(data=geneTableRemainOutliers, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor)), color='red', alpha=0.3) +  
	geom_point(data=geneTableOutliers, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor)), color='red', alpha=0.3) + 
	geom_point(data=geneTableOutliers[sampOutliers,], aes(x=log(stop), y=log(fs), colour=transcript), shape=17, size=7) +
	geom_point(data=geneTableCurated, aes(x=log(stop), y=log(fs)), color='green3', size=5) +
	theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + 
	xlab('Premature stop codons') + ylab('Frameshifts')

# fracDelBases Vs. stop codons
p2 <- ggplot() + 
	geom_point(data=geneTableRemain, aes(x=jitter(log(stop),jitterFactor), y=fracDelBases), alpha=0.3) + 
	geom_point(data=geneTable, aes(x=jitter(log(stop),jitterFactor), y=fracDelBases), color='steelblue', alpha=0.3) + 
	geom_point(data=geneTableRemainOutliers, aes(x=jitter(log(stop),jitterFactor), y=fracDelBases), color='red', alpha=0.3) +  
	geom_point(data=geneTableOutliers, aes(x=jitter(log(stop),jitterFactor), y=fracDelBases), color='red', alpha=0.3) + 
	geom_point(data=geneTableOutliers[sampOutliers,], aes(x=log(stop), y=fracDelBases, colour=transcript), shape=17, size=7) +
	geom_point(data=geneTableCurated, aes(x=log(stop), y=fracDelBases), color='green3', size=5) +
	theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + 
	xlab('Premature stop codons') + ylab('fraction of deleted bases')

# fracDelBases Vs. frameshifts
p3 <- ggplot() + 
	geom_point(data=geneTableRemain, aes(x=jitter(log(fs),jitterFactor), y=fracDelBases), alpha=0.3) + 
	geom_point(data=geneTable, aes(x=jitter(log(fs),jitterFactor), y=fracDelBases), color='steelblue', alpha=0.3) + 
	geom_point(data=geneTableRemainOutliers, aes(x=jitter(log(fs),jitterFactor), y=fracDelBases), color='red', alpha=0.3) +  
	geom_point(data=geneTableOutliers, aes(x=jitter(log(fs),jitterFactor), y=fracDelBases), color='red', alpha=0.3) + 
	geom_point(data=geneTableOutliers[sampOutliers,], aes(x=log(fs), y=fracDelBases, colour=transcript), shape=17, size=7) +
	geom_point(data=geneTableCurated, aes(x=log(fs), y=fracDelBases), color='green3', size=5) +
	theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + 
	xlab('Frameshifts') + ylab('fraction of deleted bases')

# KaKs Vs. frameshifts
p4 <- ggplot() + 
	geom_point(data=geneTableRemain, aes(x=jitter(log(fs),jitterFactor), y=KaKs), alpha=0.3) + 
	geom_point(data=geneTable, aes(x=jitter(log(fs),jitterFactor), y=KaKs), color='steelblue', alpha=0.3) + 
	geom_point(data=geneTableRemainOutliers, aes(x=jitter(log(fs),jitterFactor), y=KaKs), color='red', alpha=0.3) +  
	geom_point(data=geneTableOutliers, aes(x=jitter(log(fs),jitterFactor), y=KaKs), color='red', alpha=0.3) + 
	geom_point(data=geneTableOutliers[sampOutliers,], aes(x=log(fs), y=KaKs, colour=transcript), shape=17, size=7) +
	geom_point(data=geneTableCurated, aes(x=jitter(log(fs),jitterFactor), y=KaKs), color='green3', size=5) +
	theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + 
	xlab('Frameshifts') + ylab('Ka/Ks')

# KaKs Vs. premature stop-codons
p5 <- ggplot() + 
	geom_point(data=geneTableRemain, aes(x=jitter(log(stop),jitterFactor), y=KaKs), alpha=0.3) + 
	geom_point(data=geneTable, aes(x=jitter(log(stop),jitterFactor), y=KaKs), color='steelblue', alpha=0.3) + 
	geom_point(data=geneTableRemainOutliers, aes(x=jitter(log(stop),jitterFactor), y=KaKs), color='red', alpha=0.3) +  
	geom_point(data=geneTableOutliers, aes(x=jitter(log(stop),jitterFactor), y=KaKs), color='red', alpha=0.3) + 
	geom_point(data=geneTableOutliers[sampOutliers,], aes(x=log(stop), y=KaKs, colour=transcript), shape=17, size=7) +
	geom_point(data=geneTableCurated, aes(x=jitter(log(stop),jitterFactor), y=KaKs), color='green3', size=5) +
	theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + 
	xlab('stop') + ylab('Ka/Ks')

# KaKs Vs. premature stop-codons
p6 <- ggplot() + 
	geom_point(data=geneTableRemain, aes(x=fracDelBases, y=KaKs), alpha=0.3) + 
	geom_point(data=geneTable, aes(x=fracDelBases, y=KaKs), color='steelblue', alpha=0.3) + 
	geom_point(data=geneTableRemainOutliers, aes(x=fracDelBases, y=KaKs), color='red', alpha=0.3) +  
	geom_point(data=geneTableOutliers, aes(x=fracDelBases, y=KaKs), color='red', alpha=0.3) + 
	geom_point(data=geneTableOutliers[sampOutliers,], aes(x=fracDelBases, y=KaKs, colour=transcript), shape=17, size=7) +
	geom_point(data=geneTableCurated, aes(x=fracDelBases, y=KaKs), color='green3', size=5) +
	theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + 
	xlab('fraction of deleted bases') + ylab('Ka/Ks')

## some density plots:
d1 <- ggplot() +
	stat_density2d(data=geneTable, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor), fill=..level..), geom="polygon") + 
	scale_fill_gradient(low="lightgray", high="black", limits=c(0,3)) + 
	geom_point(data=geneTable, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor)), color='steelblue', alpha=0.3) +
	theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + 
	xlab('Premature stop codons') + ylab('Frameshifts')

d2 <- ggplot() +
	stat_density2d(data=geneTableRemain, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor), fill=..level..), geom="polygon") + 
	scale_fill_gradient(low="lightgray", high="black", limits=c(0,3)) + 
	geom_point(data=geneTableRemain, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor)), alpha=0.3) +
	theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + 
	xlab('Premature stop codons') + ylab('Frameshifts')

d3 <- ggplot() +
	stat_density2d(data=rbind(geneTableRemainOutliers, geneTableOutliers[-2]), aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor), fill=..level..), geom="polygon") + 
	scale_fill_gradient(low="lightgray", high="black", limits=c(0,3)) + 
	geom_point(data=geneTableRemainOutliers, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor)),color='red', alpha=0.3) +
	geom_point(data=geneTableOutliers, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor)),color='red', alpha=0.3) +
	theme_bw() + theme(panel.border=element_blank()) + theme(axis.text.x = element_text(vjust = 0.5, size=16, colour = "black"), axis.text.y = element_text(hjust = 1, size=16, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + 
	xlab('Premature stop codons') + ylab('Frameshifts')


pdf(paste0('output/predcitedOutliersPlots.',species,'.gamma_',gamma,'.nu_',nu,'.canonical.pdf'), onefile = TRUE)
p1
p2
p3
p4
p5
p6
d1
d2
d3
dev.off()

#ggplot() +stat_density2d(data=geneTable, aes(x=fracDelBases, y=KaKs, fill=..level..), geom="polygon")+ scale_fill_gradient(low="lightgrey", high="darkblue")

#ggplot() + geom_point(data=geneTableRemain, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor)), , alpha=0.7) + geom_point(data=geneTable, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor)), color='steelblue', alpha=0.3) + geom_point(data=geneTableCurated, aes(x=log(stop), y=log(fs)), color='green3', size=5) +  geom_point(data=geneTableRemainOutliers, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor)), color='red', alpha=0.7) +  geom_point(data=geneTableOutliers, aes(x=jitter(log(stop),jitterFactor), y=jitter(log(fs),jitterFactor)), color='red', alpha=0.7)



