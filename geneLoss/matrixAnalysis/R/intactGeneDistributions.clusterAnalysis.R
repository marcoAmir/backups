require(ggplot2)
#library("scatterplot3d")
source('/cluster/u/amirma/git/forwardGenomics/geneLoss/matrixAnalysis/src/R/multiplot.R')

species <- QUERY
geneTable <- read.table(paste0(species,'.intact.proteinCodingTx.humanOrtho'), col.names=c('assembly','qGeneID','humanGene','transcript','name','stop','fs','delExons','totalExons','fracDelExons','delCodingBases','totalCodingBases','fracDelBases','percentId','KaKs'))
geneTableRemain <- read.table(paste0(species,'.remainingCalls.proteinCodingTx.humanOrtho'), col.names=c('assembly','humanGene','transcript','name','stop','fs','delExons','totalExons','fracDelExons','delCodingBases','totalCodingBases','fracDelBases','percentId','KaKs'))


geneTableRemain[which(geneTableRemain$fracDelBases>0.9),]$KaKs <- rep(NA, length(which(geneTableRemain$fracDelBases>0.9)))
geneTable[which(geneTable$fracDelBases>0.9),]$KaKs <- rep(NA, length(which(geneTable$fracDelBases>0.9)))

curated <- read.table(paste0(species,'.curated'), col.names=c('curatedGenes'))

allGenes <- rbind(cbind(geneTable[,-2], type=rep('intact', nrow(geneTable))), cbind(geneTableRemain, type=rep('other', nrow(geneTableRemain))))
allGenes$stop             <- as.numeric(allGenes$stop)
allGenes$fs               <- as.numeric(allGenes$fs)
allGenes$delExons         <- as.numeric(allGenes$delExons)
allGenes$totalExons       <- as.numeric(allGenes$totalExons)
allGenes$delCodingBases   <- as.numeric(allGenes$delCodingBases)
allGenes$totalCodingBases <- as.numeric(allGenes$totalCodingBases)
allGenes$percentId        <- as.numeric(paste(allGenes$percentId))

allGenes$type <- as.character(allGenes$type)
allGenes[which(allGenes$humanGene %in% curated$curatedGenes), 15] <- rep('curated', length(which(allGenes$humanGene %in% curated$curatedGenes)))

curatedGenes <- subset(allGenes, humanGene %in% curated$curatedGenes)

quantileFeatures <- data.frame(percentile=0:100)	# Here I will compute the quantiles across each feature


#m <- quantile(allGenes$fs, probs=seq(0, 1, 0.001))[1000]
m <- 1.2*max(curatedGenes$fs)
p1 <- ggplot() + geom_histogram(data=geneTableRemain, aes(x=fs, y=..density..), binwidth=1, fill='red', alpha=0.5) + geom_histogram(data=geneTable, aes(x=fs, y=..density..), binwidth=1, fill='steelblue', alpha=0.5) + geom_vline(data=curatedGenes, aes(xintercept=fs, colour=name)) + scale_x_continuous(limits=c(0,m)) + theme(axis.text.x = element_text(vjust = 0.5, size=12, colour = "black"), axis.text.y = element_text(hjust = 1, size=12, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('frameshifts')

m <- 1.2*max(curatedGenes$stop)
p2 <- ggplot() + geom_histogram(data=geneTableRemain, aes(x=stop, y=..density..), binwidth=1, fill='red', alpha=0.5) + geom_histogram(data=geneTable, aes(x=stop, y=..density..), binwidth=1, fill='steelblue', alpha=0.5) + geom_vline(data=curatedGenes, aes(xintercept=stop, colour=name)) + scale_x_continuous(limits=c(0,m)) + theme(axis.text.x = element_text(vjust = 0.5, size=12, colour = "black"), axis.text.y = element_text(hjust = 1, size=12, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('stop')
quantileFeatures <- cbind(quantileFeatures, stop=quantile(allGenes$stop, probs=seq(0, 1, by=0.01)))
quantileFeatures <- cbind(quantileFeatures, fs=quantile(allGenes$fs, probs=seq(0, 1, by=0.01)))

m <- 1.2*max(curatedGenes$delExons)
p3 <- ggplot() + geom_histogram(data=geneTableRemain, aes(x=delExons, y=..density..), binwidth=1, fill='red', alpha=0.5) + geom_histogram(data=geneTable, aes(x=delExons, y=..density..), binwidth=1, fill='steelblue', alpha=0.5) + geom_vline(data=curatedGenes, aes(xintercept=delExons, colour=name)) + scale_x_continuous(limits=c(0,m)) + theme(axis.text.x = element_text(vjust = 0.5, size=12, colour = "black"), axis.text.y = element_text(hjust = 1, size=12, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('deleted-exons')
quantileFeatures <- cbind(quantileFeatures, delExons=quantile(allGenes$delExons, probs=seq(0, 1, by=0.01)))

m <- 2*max(curatedGenes$totalExons)
p4 <- ggplot() + geom_histogram(data=geneTableRemain, aes(x=totalExons, y=..density..), binwidth=1, fill='red', alpha=0.5) + geom_histogram(data=geneTable, aes(x=totalExons, y=..density..), binwidth=1, fill='steelblue', alpha=0.5) + geom_vline(data=curatedGenes, aes(xintercept=totalExons, colour=name)) + scale_x_continuous(limits=c(0,m)) + theme(axis.text.x = element_text(vjust = 0.5, size=12, colour = "black"), axis.text.y = element_text(hjust = 1, size=12, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('total-Exons')

m <- 1.1*max(curatedGenes$fracDelExons)
p5 <- ggplot() + geom_histogram(data=geneTableRemain, aes(x=fracDelExons, y=..density..), binwidth=0.05, fill='red', alpha=0.5) + geom_histogram(data=geneTable, aes(x=fracDelExons, y=..density..), binwidth=0.05, fill='steelblue', alpha=0.5) + geom_vline(data=curatedGenes, aes(xintercept=fracDelExons, colour=name)) + scale_x_continuous(limits=c(0,1)) + theme(axis.text.x = element_text(vjust = 0.5, size=12, colour = "black"), axis.text.y = element_text(hjust = 1, size=12, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('fraction of Deleted Exons')
quantileFeatures <- cbind(quantileFeatures, fracDelExons=quantile(allGenes$fracDelExons, probs=seq(0, 1, by=0.01)))

m <- 1.2*max(curatedGenes$delCodingBases)
p6 <- ggplot() + geom_histogram(data=geneTableRemain, aes(x=delCodingBases, y=..density..), binwidth=100, fill='red', alpha=0.5) + geom_histogram(data=geneTable, aes(x=delCodingBases, y=..density..), binwidth=100, fill='steelblue', alpha=0.5) + geom_vline(data=curatedGenes, aes(xintercept=delCodingBases, colour=name)) + scale_x_continuous(limits=c(0,m)) + theme(axis.text.x = element_text(vjust = 0.5, size=12, colour = "black"), axis.text.y = element_text(hjust = 1, size=12, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('deleted Coding Bases')
quantileFeatures <- cbind(quantileFeatures, delCodingBases=quantile(allGenes$delCodingBases, probs=seq(0, 1, by=0.01)))

m <- 1.2*max(curatedGenes$totalCodingBases)
p7 <- ggplot() + geom_histogram(data=geneTableRemain, aes(x=totalCodingBases, y=..density..), binwidth=100, fill='red', alpha=0.5) + geom_histogram(data=geneTable, aes(x=totalCodingBases, y=..density..), binwidth=100, fill='steelblue', alpha=0.5) + geom_vline(data=curatedGenes, aes(xintercept=totalCodingBases, colour=name)) + scale_x_continuous(limits=c(0,m)) + theme(axis.text.x = element_text(vjust = 0.5, size=12, colour = "black"), axis.text.y = element_text(hjust = 1, size=12, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('total Coding Bases')

m <- 1.2*max(curatedGenes$fracDelBases)
p8 <- ggplot() + geom_histogram(data=geneTableRemain, aes(x=fracDelBases, y=..density..), binwidth=0.05, fill='red', alpha=0.5) + geom_histogram(data=geneTable, aes(x=fracDelBases, y=..density..), binwidth=0.05, fill='steelblue', alpha=0.5) + geom_vline(data=curatedGenes, aes(xintercept=fracDelBases, colour=name)) + scale_x_continuous(limits=c(0,1)) + theme(axis.text.x = element_text(vjust = 0.5, size=12, colour = "black"), axis.text.y = element_text(hjust = 1, size=12, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('fraction of Deleted Bases')
quantileFeatures <- cbind(quantileFeatures, fracDelBases=quantile(allGenes$fracDelBases, probs=seq(0, 1, by=0.01)))

curatedGenes$percentId <- as.numeric(paste(curatedGenes$percentId))
m <- 2*max(as.numeric(paste(curatedGenes$percentId)), na.rm=TRUE)
p9 <- ggplot() + geom_histogram(data=geneTableRemain[which(geneTableRemain$percentId!="?"),], aes(x=as.numeric(paste(percentId)), y=..density..), binwidth=0.05, fill='red', alpha=0.5) + geom_histogram(data=geneTable[which(geneTable$percentId!="?"),], aes(x=as.numeric(paste(percentId)), y=..density..), binwidth=0.05, fill='steelblue', alpha=0.5) + geom_vline(data=curatedGenes, aes(xintercept=percentId, colour=name)) + scale_x_continuous(limits=c(0,m)) + theme(axis.text.x = element_text(vjust = 0.5, size=12, colour = "black"), axis.text.y = element_text(hjust = 1, size=12, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('percent Id')

m <- 2*max(curatedGenes$KaKs, na.rm=TRUE)
p10 <- ggplot() + geom_histogram(data=geneTableRemain, aes(x=KaKs, y=..density..), binwidth=0.01, fill='red', alpha=0.5) + geom_histogram(data=geneTable, aes(x=KaKs, y=..density..), binwidth=0.01, fill='steelblue', alpha=0.5) + geom_vline(data=curatedGenes, aes(xintercept=KaKs, colour=name)) + scale_x_continuous(limits=c(0,m)) + theme(axis.text.x = element_text(vjust = 0.5, size=12, colour = "black"), axis.text.y = element_text(hjust = 1, size=12, colour = "black"), axis.title=element_text(size=12),axis.title.y=element_text(vjust=0.3)) + xlab('Ka/Ks')
quantileFeatures <- cbind(quantileFeatures, KaKs=quantile(allGenes$KaKs, probs=seq(0, 1, by=0.01), na.rm = TRUE))

write.table(quantileFeatures, file = paste0(species, ".Quantiles"), sep="\t", quote = FALSE, row.names = FALSE)

pdf(paste0('oneDimensionalDists.',species,'.pdf'), onefile = TRUE)
multiplot(p1, p3, p4, p8, p9, p2, p6, p7, p5, p10, cols=2)

# PCA
circle <- function(center = c(0, 0), npoints = 100) {
    r = 1
    tt = seq(0, 2 * pi, length = npoints)
    xx = center[1] + r * cos(tt)
    yy = center[1] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}

forPCA <- allGenes[,5:14]
rownames(forPCA) <- allGenes$humanGene
pca <- prcomp(~stop+fs+delExons+totalExons+fracDelExons+delCodingBases+totalCodingBases+fracDelBases+percentId+KaKs, data=forPCA, center= TRUE, scale = TRUE, na.action = na.omit)
allGenesPCA <- as.data.frame(pca$x)
allGenesPCA <- cbind(humanGene = rownames(allGenesPCA), allGenesPCA)
allGenesPCA <- merge(x = allGenesPCA, y = allGenes, by = 'humanGene')[,c(1:10, 14, 25, 12)]

p11 <- ggplot() + geom_hline(yintercept = 0, colour="gray65") + geom_vline(xintercept = 0, colour="gray65") + geom_point(data = subset(allGenesPCA, type=="intact"), aes(x=PC1, y=PC2, colour=type, alpha=0.6)) + geom_point(data = subset(allGenesPCA, type=="other"), aes(x=PC1, y=PC2, colour=type, alpha=0.3)) + geom_point(data = subset(allGenesPCA, type=="curated"), aes(x=PC1, y=PC2, colour=type), size=4) + scale_color_manual(values=c('green3','steelblue','red')) + ggtitle(paste0(species," - PCA of gene loss features"))

p12 <- ggplot() + geom_hline(yintercept = 0, colour="gray65") + geom_vline(xintercept = 0, colour="gray65") + geom_point(data = subset(allGenesPCA, type=="intact"), aes(x=PC1, y=PC3, colour=type, alpha=0.6)) + geom_point(data = subset(allGenesPCA, type=="other"), aes(x=PC1, y=PC3, colour=type, alpha=0.3)) + geom_point(data = subset(allGenesPCA, type=="curated"), aes(x=PC1, y=PC3, colour=type), size=4) + scale_color_manual(values=c('green3','steelblue','red'))

multiplot(p11, p12, cols=1)


# post-pca:
axisRange <- function(d){
	u <- seq(from=min(d), to=max(d), by=(max(d)-min(d))/20)
	return(u)
	
}

getFraction <- function(d, speciesData){
	axisDens <- data.frame()
	for(t in d){
		a <- speciesData[which(speciesData[,1]>t), 2]
		other <- length(which(a=="other"))
		curated <- length(which(a=="curated"))
		intact <- length(which(a=="intact"))
		axisDens <- rbind(axisDens, data.frame(intact = intact/sum(other, curated, intact), other = other/sum(other, curated, intact), curated = curated/sum(other, curated, intact)))
		#axisDens <- rbind(axisDens, data.frame(intact = intact, other = other, curated = curated))
	}
	return(axisDens)
}

ax <- data.frame(apply(allGenesPCA[,2:5], 2, FUN = axisRange))
pcAxisDens <- data.frame()
for (c in 1:ncol(ax)){
	axisDens <- getFraction(ax[,c], allGenesPCA[,c(c+1, 12)])
	pcAxisDens <- rbind(pcAxisDens, cbind(axis = ax[,c], axisDens, principalComponent=rep(colnames(ax)[c])))
}
ggplot(pcAxisDens) + geom_line(aes(x = axis, y= intact), color="steelblue") + geom_line(aes(x = axis, y= other), color="red") + geom_line(aes(x = axis, y= curated), color="green3") + facet_wrap(~principalComponent, scales = "free")

# Outlier analysis:
getMeanPercentile <- function(featureRead, c){
	if(sum(quantileFeatures[, c+1] >= featureRead)<101 & sum(quantileFeatures[, c+1] >= featureRead)>0){
		return(which((1:101)*(quantileFeatures[, c+1] >= featureRead)!=0)[1] - 1)
	} else {
		return((sum(quantileFeatures[, c+1] >= featureRead)==101)*0 + (sum(quantileFeatures[, c+1] >= featureRead)==0)*100)
	}	
}

outlier <- c()

for (i in 1:nrow(allGenes)){
	featurePercentile <- c()
	c <- 1
	for (j in c(5:7, 9:10, 12, 14)){
		if(!is.na(allGenes[i,j])){
			featurePercentile <- c(featurePercentile, getMeanPercentile(allGenes[i,j], c))
			c <- c + 1
		}
	}
	outlier <- c(outlier,  mean(featurePercentile))
}

allGenes <- cbind(allGenes, outlier=outlier)

ggplot() + geom_violin(data=subset(allGenes, type=="intact"), aes(x=1, y=outlier), fill="steelblue") + geom_violin(data=subset(allGenes, type=="other"), aes(x=2, y=outlier), fill="red") + geom_hline(data=subset(allGenes, type=="curated"), aes(yintercept=outlier, colour=name)) + geom_text(data=subset(allGenes, type=="curated"), aes(x=jitter(rep(1.5, nrow(subset(allGenes, type=="curated"))), 30), y=subset(allGenes, type=="curated")$outlier), label=subset(allGenes, type=="curated")$name)

dev.off()
