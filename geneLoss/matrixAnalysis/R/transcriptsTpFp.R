require(ggplot2)

tpGenes <- read.table('TP/tpGeneTranscripts.tsv', header=TRUE)
fpGenes <- read.table('FP/fpGeneTranscripts.tsv', header=TRUE)

ggplot() + geom_point(data=tpGenes, aes(x=num.total.exons, y=num.deleted.exons, color='TP'), alpha=0.5, size=5) + geom_point(data=fpGenes, aes(x=num.total.exons, y=num.deleted.exons, color='FP'), alpha=0.5, size=5) + xlab('Total transcript exons') + ylab('Deleted transcript exons') + ggtitle('currSpecies, deleted exons') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "blue"), labels=c("FP", "TP"))
ggsave('exonsDeletions.png')

ggplot() + geom_point(data=tpGenes, aes(x=num.total.bases, y=num.bases.in.deleted.exons, color='TP'), alpha=0.5, size=5) + geom_point(data=fpGenes, aes(x=num.total.bases, y=num.bases.in.deleted.exons, color='FP'), alpha=0.5, size=5) + xlab('Total bases in transcript exons') + ylab('Deleted bases in transcript exons') + ggtitle('currSpecies, deleted bases') + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=20), axis.title=element_text(size=20)) + theme(panel.border = element_blank()) + scale_colour_manual(name='gene loss call', values=c("red", "blue"), labels=c("FP", "TP"))
ggsave('exonsBasesDeletions.png')

ggplot() + geom_point(data=tpGenes, aes(x=1, y=fraction.deleted.exons, color='TP'), alpha=0.5, size=5) + geom_point(data=fpGenes, aes(x=1, y=fraction.deleted.exons, color='FP'), alpha=0.5, size=5) + geom_point(data=tpGenes, aes(x=2, y=fraction.deleted.bases, color='TP'), alpha=0.5, size=5) + geom_point(data=fpGenes, aes(x=2, y=fraction.deleted.bases, color='FP'), alpha=0.5, size=5) + theme_bw() + scale_x_continuous(limits=c(0.5,2.5), breaks=c(1,2), labels=c('Frac. of del. exons','Frac. of del. bases')) + scale_colour_manual(name='gene loss call', values=c("red", "blue"), labels=c("FP", "TP")) + xlab('') + ylab('fraction') + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=14), axis.title=element_text(size=14)) + theme(panel.border = element_blank()) + ggtitle('currSpecies, del. exons and del. bases fractions')
ggsave('fractionDeletions.png')


# exonsDel <- data.frame(ensgID=tpGenes$ensgID, num.deleted.exons=tpGenes$num.deleted.exons, num.total.exons=tpGenes$num.total.exons, call=rep('TP', nrow(tpGenes)))


