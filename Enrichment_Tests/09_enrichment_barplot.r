#PURPOSE: Generate barplot comparing fold enrichment of GO category genes for different datasets
#~/software/envs/r4.3.3/ conda environment has ggpattern installed

library(ggplot2)
library(ggpattern)

#Read in GO enrichment results from continuous RERconverge analyses
#OLD - based on 250 most decelerated
#continuous_eutheria<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/RERconverge.continuous.fullEuth.GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.csv", header=TRUE)
#continuous_laurasiatheria<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/RERconverge.continuous.laurasubset.GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.csv", header=TRUE)
#continuous_rodent<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/RERconverge.continuous.rodentsubset.GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.csv", header=TRUE)
#continuous_primate<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/RERconverge.continuous.primatesubset.GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.csv", header=TRUE)
#NEW - based on all significant (raw P < 0.05) genes
continuous_eutheria<-read.csv("./GOresults_Sig/RERconverge_GOenrichment.sigdecOfAll.translated.keepEarlyStops.noAsk.csv", header=TRUE)
continuous_laurasiatheria<-read.csv("./GOresults_Sig/RERconverge.LAURA_GOenrichment.sigdecOfAll.translated.keepEarlyStops.noAsk.csv", header=TRUE)
continuous_rodent<-read.csv("./GOresults_Sig/RERconverge.RODENT_GOenrichment.sigdecOfAll.translated.keepEarlyStops.noAsk.csv", header=TRUE)
continuous_primate<-read.csv("./GOresults_Sig/RERconverge.PRIMATE_GOenrichment.sigdecOfAll.translated.keepEarlyStops.noAsk.csv", header=TRUE)

#Read in GO enrichment results from binary RERconverge analyses
#DECELERATED - WRONG
#binary_eutheria<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/binary_enrich/250_genes/RERconverge.LESS.STRICT.binary.noCOMP._GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.noAsk.csv", header=TRUE)
#binary_laurasiatheria<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/binary_enrich/250_genes/RERconverge.LAURA.binary.noCOMP._GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.noAsk.csv", header=TRUE)
#binary_rodent<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/binary_enrich/250_genes/RERconverge.RODENT.binary.noCOMP._GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.noAsk.csv", header=TRUE)
#binary_primate<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/binary_enrich/250_genes/RERconverge.PRIMATE.binary.noCOMP._GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.noAsk.csv", header=TRUE)
#ACCELERATED - CORRECT
#binary_eutheria<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/binary_enrich/250_genes/RERconverge.LESS.STRICT.binary.noCOMP._GOenrichment.top250accOfAll.translated.keepEarlyStops.NoTraitWinsor.noAsk.csv", header=TRUE)
#binary_laurasiatheria<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/binary_enrich/250_genes/RERconverge.LAURA.binary.noCOMP._GOenrichment.top250accOfAll.translated.keepEarlyStops.NoTraitWinsor.noAsk.csv", header=TRUE)
#binary_rodent<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/binary_enrich/250_genes/RERconverge.RODENT.binary.noCOMP._GOenrichment.top250accOfAll.translated.keepEarlyStops.NoTraitWinsor.noAsk.csv", header=TRUE)
#binary_primate<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/binary_enrich/250_genes/RERconverge.PRIMATE.binary.noCOMP._GOenrichment.top250accOfAll.translated.keepEarlyStops.NoTraitWinsor.noAsk.csv", header=TRUE)
#NEW - based on all significant (raw P < 0.05) genes
binary_eutheria<-read.csv("./GOresults_Sig/RERconverge.LESS.STRICT.binary.noCOMP._GOenrichment.sigaccOfAll.translated.keepEarlyStops.noAsk.csv", header=TRUE)
binary_laurasiatheria<-read.csv("./GOresults_Sig/RERconverge.LAURA.binary.noCOMP._GOenrichment.sigaccOfAll.translated.keepEarlyStops.noAsk.csv", header=TRUE)
binary_rodent<-read.csv("./GOresults_Sig/RERconverge.RODENT.binary.noCOMP._GOenrichment.sigaccOfAll.translated.keepEarlyStops.noAsk.csv", header=TRUE)
binary_primate<-read.csv("./GOresults_Sig/RERconverge.PRIMATE.binary.noCOMP._GOenrichment.sigaccOfAll.translated.keepEarlyStops.noAsk.csv", header=TRUE)

#List of datasets
datasets<-c("continuous_eutheria", "continuous_laurasiatheria", "continuous_rodent", "continuous_primate", "binary_eutheria", "binary_laurasiatheria", "binary_rodent", "binary_primate")

#List of GO terms we want to plot
go_terms<-c("spermatogenesis", "sperm motility", "oogenesis", "gamete generation", "meiotic cell cycle", "sexual reproduction", "cilium movement", "spermatid differentiation", "fertilization", "sensory perception of smell") #"reproductive process",

#Make a data frame for ggplot2 input
#GOterm	fold_enrich	dataset	analysis_type
mydf<-c()
for(i in go_terms){
	for(j in datasets){
		ds<-get(j)
		sig<-ds$Significant[which(ds$Term==i)]
		exp<-ds$Expected[which(ds$Term==i)]
		newline<-c(i, sig/exp, sub(".*_","",j), sub("_.*","",j))
		#print(newline)
		mydf<-rbind(mydf, newline)
	}
}

mydf<-as.data.frame(mydf)
colnames(mydf)<-c("GOterm", "Fold_Enrichment", "Dataset", "Analysis_Type")
mydf$Fold_Enrichment<-as.numeric(as.character(mydf$Fold_Enrichment))

#Plot
p<-ggplot(data=mydf, aes(x=GOterm, y=Fold_Enrichment, fill=Dataset, pattern=Analysis_Type))
p<-p + geom_bar_pattern(stat="identity", position="dodge", pattern_spacing = 0.01, pattern_frequency = 5, pattern_angle = 45)
p<-p + geom_hline(yintercept=1) + theme_minimal()
p<-p + scale_fill_manual(values=c("#56B4E9", "#009E73", "#E69F00", "#CC79A7"))
p<-p + scale_pattern_manual(values=c("none", "stripe"))

pdf("enrichment_barplots.sigDeceleratedGenes.accInBinary.pdf", height=7, width=14)
print(p)
dev.off()

#Make a smaller plot for publication with fewer GO terms
#smalldf<-mydf[which(mydf$GOterm %in% c("sperm motility", "spermatid differentiation", "spermatogenesis", "sensory perception of smell")),]
smalldf<-mydf[which(mydf$GOterm %in% c("sperm motility", "spermatid differentiation", "spermatogenesis")),]
p<-ggplot(data=smalldf, aes(x=GOterm, y=Fold_Enrichment, fill=Dataset, pattern=Analysis_Type))
p<-p + geom_bar_pattern(stat="identity", position="dodge", pattern_spacing = 0.01, pattern_frequency = 5, pattern_angle = 45, color="black")
p<-p + geom_hline(yintercept=1) + theme_minimal()
p<-p + scale_fill_manual(values=c("#56B4E9", "#009E73", "#E69F00", "#CC79A7"))
p<-p + scale_pattern_manual(values=c("none", "stripe"))

pdf("enrichment_barplots.sigDeceleratedGenes.accInBinary.FORPUBLICATION.pdf")
print(p)
dev.off()

#Make smaller plot that includes smell term
smalldf<-mydf[which(mydf$GOterm %in% c("sperm motility", "sensory perception of smell", "spermatogenesis")),]
p<-ggplot(data=smalldf, aes(x=GOterm, y=Fold_Enrichment, fill=Dataset, pattern=Analysis_Type))
p<-p + geom_bar_pattern(stat="identity", position="dodge", pattern_spacing = 0.01, pattern_frequency = 5, pattern_angle = 45, color="black")
p<-p + geom_hline(yintercept=1) + theme_minimal()
p<-p + scale_fill_manual(values=c("#56B4E9", "#009E73", "#E69F00", "#CC79A7"))
p<-p + scale_pattern_manual(values=c("none", "stripe"))

pdf("enrichment_barplots.sigDeceleratedGenes.accInBinary.WITHSMELL.pdf")
print(p)
dev.off()

