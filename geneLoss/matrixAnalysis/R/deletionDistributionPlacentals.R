require(ggplot2)

speciesTable <- data.frame(species=c("Chimp","Gorilla","Orangutan","Gibbon","Rhesus","Marmoset","Bushbaby","Tarsier"), group=rep("Primates", 8))				# Add the primates
speciesTable <- rbind(speciesTable, data.frame(species=c("Treeshrew","Squirrel","Mouse","Rat","KangarooRat","GuineaPig","Rabbit","Pika"), group=rep("Rodents", 8)))		# Add the rodents
speciesTable <- rbind(speciesTable, data.frame(species=c("Alpaca","Dolphin","Cow","Sheep","Horse","Dog","Cat","Microbat","Megabat","Shrew","Hedgehog","Panda","Elephant","Tenrec"), group=rep("Placentals", 14)))		# add other placentals (Laura+Afro)
#speciesTable <- rbind(speciesTable, data.frame(species=c("Elephant","Tenrec"), group=rep("Afrotheria", 2)))		# add other afrotherians

tpGenes <- data.frame()
fpGenes <- data.frame()

for(i in 1:nrow(speciesTable)){
	currGroup <- as.character(speciesTable$group[i])
	tpTmp <- read.table(paste(as.character(speciesTable$species[i]),'/TP/tpGeneTranscripts.tsv', sep=""), header=TRUE)
	fpTmp <- read.table(paste(as.character(speciesTable$species[i]),'/FP/fpGeneTranscripts.tsv', sep=""), header=TRUE)
	tpGenes   <- rbind(tpGenes, cbind(tpTmp, group=rep(currGroup, nrow(tpTmp))))
	fpGenes   <- rbind(fpGenes, cbind(fpTmp, group=rep(currGroup, nrow(fpTmp))))
}

##### PLOTS: #######

dir.create(paste(getwd(),'/featuresTpFpCalls', sep=""), recursive=TRUE, showWarnings=FALSE)

#### 1. fraction of deleted exons
# primates:
ggplot() + geom_density(data=subset(tpGenes, group=='Primates'), aes(x=fraction.deleted.exons, color='TP')) + geom_density(data=subset(fpGenes, group=='Primates'), aes(x=fraction.deleted.exons, color='FP')) + xlab('Fraction of deleted exons') + ggtitle('Primates, TP/FP deleted exons fraction') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "darkblue"), labels=c("FP", "TP")) + scale_x_continuous(limits=c(0,1))
ggsave('featuresTpFpCalls/primates.delExonsFractions.png')
# rodents:
ggplot() + geom_density(data=subset(tpGenes, group=='Rodents'), aes(x=fraction.deleted.exons, color='TP')) + geom_density(data=subset(fpGenes, group=='Rodents'), aes(x=fraction.deleted.exons, color='FP')) + xlab('Fraction of deleted exons') + ggtitle('Rodents, TP/FP deleted exons fraction') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "darkblue"), labels=c("FP", "TP")) + scale_x_continuous(limits=c(0,1))
ggsave('featuresTpFpCalls/rodents.delExonsFractions.png')
# placentals/afro:
ggplot() + geom_density(data=subset(tpGenes, group=='Placentals'), aes(x=fraction.deleted.exons, color='TP')) + geom_density(data=subset(fpGenes, group=='Placentals'), aes(x=fraction.deleted.exons, color='FP')) + xlab('Fraction of deleted exons') + ggtitle('Placentals, TP/FP deleted exons fraction') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "darkblue"), labels=c("FP", "TP")) + scale_x_continuous(limits=c(0,1))
ggsave('featuresTpFpCalls/placentals.delExonsFractions.png')
# All the species:
ggplot() + geom_density(data=tpGenes, aes(x=fraction.deleted.exons, color='TP')) + geom_density(data=fpGenes, aes(x=fraction.deleted.exons, color='FP')) + xlab('Fraction of deleted exons') + ggtitle('TP/FP deleted exons fraction') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "darkblue"), labels=c("FP", "TP")) + scale_x_continuous(limits=c(0,1))
ggsave('featuresTpFpCalls/all.delExonsFractions.png')

#### 2. fraction of deleted bases
# primates:
ggplot() + geom_density(data=subset(tpGenes, group=='Primates'), aes(x=fraction.deleted.bases, color='TP')) + geom_density(data=subset(fpGenes, group=='Primates'), aes(x=fraction.deleted.bases, color='FP')) + xlab('Fraction of deleted bases') + ggtitle('Primates, TP/FP deleted bases fraction') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "darkblue"), labels=c("FP", "TP")) + scale_x_continuous(limits=c(0,1))
ggsave('featuresTpFpCalls/primates.delBasesFractions.png')
# rodents:
ggplot() + geom_density(data=subset(tpGenes, group=='Rodents'), aes(x=fraction.deleted.bases, color='TP')) + geom_density(data=subset(fpGenes, group=='Rodents'), aes(x=fraction.deleted.bases, color='FP')) + xlab('Fraction of deleted bases') + ggtitle('Rodents, TP/FP deleted bases fraction') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "darkblue"), labels=c("FP", "TP")) + scale_x_continuous(limits=c(0,1))
ggsave('featuresTpFpCalls/rodents.delBasesFractions.png')
# placentals/afro:
ggplot() + geom_density(data=subset(tpGenes, group=='Placentals'), aes(x=fraction.deleted.bases, color='TP')) + geom_density(data=subset(fpGenes, group=='Placentals'), aes(x=fraction.deleted.bases, color='FP')) + xlab('Fraction of deleted bases') + ggtitle('Placentals, TP/FP deleted bases fraction') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "darkblue"), labels=c("FP", "TP")) + scale_x_continuous(limits=c(0,1))
ggsave('featuresTpFpCalls/placentals.delBasesFractions.png')
# All the species:
ggplot() + geom_density(data=tpGenes, aes(x=fraction.deleted.bases, color='TP')) + geom_density(data=fpGenes, aes(x=fraction.deleted.bases, color='FP')) + xlab('Fraction of deleted bases') + ggtitle('TP/FP deleted bases fraction') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "darkblue"), labels=c("FP", "TP")) + scale_x_continuous(limits=c(0,1))
ggsave('featuresTpFpCalls/all.delBasesFractions.png')

#### 3. number of deleted exons
# primates:
ggplot() + geom_density(data=subset(tpGenes, group=='Primates'), aes(x=num.deleted.exons, color='TP')) + geom_density(data=subset(fpGenes, group=='Primates'), aes(x=num.deleted.exons, color='FP')) + xlab('# deleted exons') + ggtitle('Primates, TP/FP deleted exons') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "darkblue"), labels=c("FP", "TP")) 
ggsave('featuresTpFpCalls/primates.delExons.png')
# rodents:
ggplot() + geom_density(data=subset(tpGenes, group=='Rodents'), aes(x=num.deleted.exons, color='TP')) + geom_density(data=subset(fpGenes, group=='Rodents'), aes(x=num.deleted.exons, color='FP')) + xlab('# deleted exons') + ggtitle('Rodents, TP/FP deleted exons') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "darkblue"), labels=c("FP", "TP")) 
ggsave('featuresTpFpCalls/rodents.delExons.png')
# placentals/afro:
ggplot() + geom_density(data=subset(tpGenes, group=='Placentals'), aes(x=num.deleted.exons, color='TP')) + geom_density(data=subset(fpGenes, group=='Placentals'), aes(x=num.deleted.exons, color='FP')) + xlab('# deleted exons') + ggtitle('Placentals, TP/FP deleted exons') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "darkblue"), labels=c("FP", "TP")) 
ggsave('featuresTpFpCalls/placentals.delExons.png')
# All the species:
ggplot() + geom_density(data=tpGenes, aes(x=num.deleted.exons, color='TP')) + geom_density(data=fpGenes, aes(x=num.deleted.exons, color='FP')) + xlab('# deleted exons') + ggtitle('TP/FP deleted exons') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "darkblue"), labels=c("FP", "TP")) 
ggsave('featuresTpFpCalls/all.delExons.png')
