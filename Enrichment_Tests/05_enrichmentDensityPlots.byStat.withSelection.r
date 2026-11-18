#PURPOSE: Generate plots of enrichment density by RERconverge rank

library(topGO)
#library(biomaRt)
library(ggplot2)
source("MariaBarcodePlotFinal.colorBySelec.R")

signedPMat<-c()

res0<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/RTS_RER_translated_keepEarlyStops_noAsk_results.csv", header=TRUE)
#Remove NAs
res<-res0[which(!(is.na(res0$Rho))),]
rhosign<-sign(res$Rho)
res$signLogP<-(-rhosign*log10(res$P))

signedPMat<-cbind(signedPMat,res$signLogP)
rownames(signedPMat)<-res$X

print("signedPMat info:")
colnames(signedPMat)<-c("241cont")
print(paste("Matrix?", is.matrix(signedPMat)))
print(dim(signedPMat))
print(head(signedPMat))

print("Getting gene sets...")
#Load GO data from previous topGO runs
load("GOdata_ensemblHs.rds")
GO_spermatogenesis<-unlist(genesInTerm(myGOdata, "GO:0007283"))
humanAtlasTestisSpecific<-read.table("tissue_category_rna_testis_DetectedOnlyInTestis.tsv", sep="\t", header=TRUE)
GO_spermatogenesis_ts<-intersect(GO_spermatogenesis, humanAtlasTestisSpecific$Ensembl)

print("Getting selection results...")
relax_table<-read.csv("HYPHY_RUNS/relax_earlystopfix_noAsk_results.combinedGOspermatogenesisGenes.csv", header=TRUE)
#Correct for multiple tests
relax_table$padj<-p.adjust(relax_table$p.value, method="BH")
#relax<-relax_table[intersect(which(relax_table$padj < 0.05), which(relax_table$K < 1)),]
relax<-relax_table[intersect(which(relax_table$p.value < 0.05), which(relax_table$K < 1)),]
relax_genes<-unique(sort(intersect(relax$id, rownames(signedPMat))))
print(relax_genes)
#Busted-ph P-values already corrected for multiple tests; see BUSTED-PH preprint page 4 under 'Hypothesis testing and implementation'
bustedph<-scan("HYPHY_RUNS/busted-ph_sigGenes.strict.txt", what=character())
bustedph_genes<-unique(sort(intersect(bustedph, rownames(signedPMat))))
selec_status<-c()
#1 = relaxed purifying selection; 2 = branch-site positive selection; 0 = none
for(i in rownames(signedPMat)){
        if(i %in% relax_genes){
                selec_status<-c(selec_status, 1)
        } else if(i %in% bustedph_genes){
                selec_status<-c(selec_status, 2)
        } else{
                selec_status<-c(selec_status, 0)
        }
}
signedPMat_selec<-as.matrix(cbind(signedPMat, selec_status))
print(head(signedPMat_selec))

print("Plotting with selection results...")
myCol<-c("#000000", "#CCCCCC", "#028090", "#984447")
#none = #CCCCCC, busted-ph = #984447, relaxed = #028090, enrichment line = #000000
pdf("enrichment_density_withSelection.noBHcorrection.241continuous.translated.span0_2.byStat.keepEarlyStops.noAsk.pdf", height=5, width=7, onefile=TRUE)
makeMultiGSEAplotSelection(signedPMat_selec, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis", cols=myCol)

print("Using genes with both selection information and in spermatogenesis GO category as input genes")
GO_spermatogenesis_bustedph<-intersect(rownames(signedPMat_selec)[which(signedPMat_selec[,'selec_status']==2)], GO_spermatogenesis)
print(length(GO_spermatogenesis_bustedph))
makeMultiGSEAplotSelection(signedPMat_selec, GO_spermatogenesis_bustedph, span=.2, title="GO category: spermatogenesis - selection associated with small testes", cols=c("#984447","#984447"))
GO_spermatogenesis_relaxed<-intersect(rownames(signedPMat_selec)[which(signedPMat_selec[,'selec_status']==1)], GO_spermatogenesis)
print(length(GO_spermatogenesis_relaxed))
makeMultiGSEAplotSelection(signedPMat_selec, GO_spermatogenesis_relaxed, span=.2, title="GO category: spermatogenesis - relaxed selection", cols=c("#028090","#028090"))

print("Testis-specific only")
makeMultiGSEAplotSelection(signedPMat_selec, GO_spermatogenesis_ts, span=.2, title="GO category: spermatogenesis AND testis-specific", cols=myCol)

dev.off()

