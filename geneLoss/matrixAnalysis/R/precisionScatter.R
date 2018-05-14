require(ggplot2)

hdefaultFiltering <- 'h0'	# I count valid ensembl protein coding orthologs using various filters on the human and ortholog species
odefaultFiltering <- 'o0'	# see:       /cluster/u/amirma/geneLoss/matrixAnalysis/humanFilterTypes.txt
				#       and  /cluster/u/amirma/geneLoss/matrixAnalysis/placentalMammalLosses/orthoFilterTypes.txt
				# I'll generate two scatter plots:
				#	 I - with the selected default filter on the human and orthocolored by species phylo group
				#	II - with the all the sets colored by the combinations of filters used

calls <- read.table('callsVsEnsembl.tsv', header=TRUE)
#ggplot() + geom_point(data=calls, aes(x=lossCalls, y=ensemblProtCodOrtho)) + geom_text(data=calls, aes(x=lossCalls, y=ensemblProtCodOrtho, label=species), vjust=1.1) + geom_line(data=data.frame(X=c(0,100), Y=c(0,100)), aes(x=X, y=Y), linetype=2) + theme_bw() + xlab('Gene loss calls') + ylab('Ensembl protein-coding ortho.')

x <- calls$lossCalls
y <- calls$ensemblProtCodOrtho
lmPrec <- lm(formula = y ~ x)
a <- lmPrec$coefficients[2]
b <- lmPrec$coefficients[1]
u <- 1:100
v <- a*u+b

ggplot() + geom_point(data=subset(calls, humanFilters==hdefaultFiltering & orthoFilters==odefaultFiltering), aes(x=lossCalls, y=ensemblProtCodOrtho)) + geom_text(data=subset(calls, humanFilters==hdefaultFiltering & orthoFilters==odefaultFiltering), aes(x=lossCalls, y=ensemblProtCodOrtho, label=species, color=group), vjust=1.1) + geom_line(data=data.frame(X=c(0,100), Y=c(0,100)), aes(x=X, y=Y), linetype=2) + geom_line(data=data.frame(X=u, Y=v), aes(x=X, y=Y), linetype=2) + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=24), axis.title=element_text(size=24)) + xlab('Gene loss calls') + ylab('Ensembl protein-coding ortho.')
ggsave('callsVsEnsembl.png', height=8, width=12)

colorSet <- c("darkgrey", "darkgreen", "skyblue", "red", "yellow","cyan")
ggplot() + geom_point(data=calls, aes(x=lossCalls, y=ensemblProtCodOrtho)) + geom_text(data=calls, aes(x=lossCalls, y=ensemblProtCodOrtho, label=species, color=orthoFilters), vjust=1.1) + geom_line(data=data.frame(X=c(0,100), Y=c(0,100)), aes(x=X, y=Y), linetype=2) + geom_line(data=data.frame(X=u, Y=v), aes(x=X, y=Y), linetype=2) + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=24), axis.title=element_text(size=24)) + xlab('Gene loss calls') + ylab('Ensembl protein-coding ortho.') + scale_colour_manual(values=colorSet[1:length(unique(calls$orthoFilters))])
ggsave('callsVsEnsembl.byOrthoFilters.png', height=8, width=12)

ggplot() + geom_point(data=calls, aes(x=lossCalls, y=ensemblProtCodOrtho)) + geom_text(data=calls, aes(x=lossCalls, y=ensemblProtCodOrtho, label=species, color=group), vjust=1.1) + geom_line(data=data.frame(X=c(0,100), Y=c(0,100)), aes(x=X, y=Y), linetype=2)+ geom_line(data=data.frame(X=u, Y=v), aes(x=X, y=Y), linetype=2) + theme_bw() + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.text=element_text(size=24), axis.title=element_text(size=24)) + xlab('Gene loss calls') + ylab('Ensembl protein-coding ortho.') + facet_grid(humanFilters~orthoFilters)
ggsave('callsVsEnsembl.GridHumanOrthoFilters.png', height=8, width=16)
