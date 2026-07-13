#PURPOSE: Look for correlations and overlaps among RERconverge results for all eutherians vs primates vs rodents vs laurasiatheria

library(VennDiagram)
library(topGO)

#Read in data
#eut<-read.csv("translated241_results_log_rts_csv_trusted_241continuous.csv", header=TRUE)
#pri<-read.csv("241CONTINUOUS_PRIMATES_ONLY/translated_43species_breadTree_LOGGED_MAREL_primate.csv", header=TRUE)
#rod<-read.csv("241CONTINUOUS_RODENTS_ONLY/translated_breadTree_53RODENT_RERresults.csv", header=TRUE)
#lau<-read.csv("241CONTINUOUS_LAURASIATHERIA_ONLY/translated_breadTree_lauras_only_RERresults.csv", header=TRUE)
eut<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/RTS_RER_translated_keepEarlyStops_noAsk_results.csv", header=TRUE)
pri<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/translated_breadTree_primates_only_KeepEarlyStops_noAsk_RERresults.csv", header=TRUE)
rod<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/translated_breadTree_rodents_only_KeepEarlyStops_noAsk_RERresults.csv", header=TRUE)
lau<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/translated_breadTree_lauras_only_KeepEarlyStops_noAsk_RERresults.csv", header=TRUE)

#Get significantly decelerated genes
eut_sig_dec<-eut[intersect(which(eut$Rho < 0), which(eut$P < 0.05)),]
pri_sig_dec<-pri[intersect(which(pri$Rho < 0), which(pri$P < 0.05)),]
rod_sig_dec<-rod[intersect(which(rod$Rho < 0), which(rod$P < 0.05)),]
lau_sig_dec<-lau[intersect(which(lau$Rho < 0), which(lau$P < 0.05)),]

#Generate Venn Diagram of significantly decelerated genes
print("Making Venn Diagram...")
catnames<-c("Eutherian", "Primate", "Rodent", "Laurasiatherian")
#vd<-venn.diagram(x=list(eut_sig_dec$X, pri_sig_dec$X, rod_sig_dec$X, lau_sig_dec$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.allSigDec.all241datasets.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")
vd<-venn.diagram(x=list(eut_sig_dec$X, pri_sig_dec$X, rod_sig_dec$X, lau_sig_dec$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.allSigDec.all241datasets.keepEarlyStops.noTraitWinsor.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")

catnames<-c("Primate", "Rodent", "Laurasiatherian")
#vd<-venn.diagram(x=list(eut_sig_dec$X, pri_sig_dec$X, rod_sig_dec$X, lau_sig_dec$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.allSigDec.all241datasets.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")
vd<-venn.diagram(x=list(pri_sig_dec$X, rod_sig_dec$X, lau_sig_dec$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5), filename="./enrich_vis/venn_diagram.allSigDec.all241datasetsNO_EUTH.keepEarlyStops.noTraitWinsor.png", imagetype="png", output=TRUE, fill=c("#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")



print("Genes significant in both eutherians and primates:")
print(intersect(eut_sig_dec$X, pri_sig_dec$X))
print("Genes significant in both eutherians and rodents:")
print(intersect(eut_sig_dec$X, rod_sig_dec$X))
print("Genes significant in both eutherians and laurasiatherians:")
print(intersect(eut_sig_dec$X, lau_sig_dec$X))
print("Genes significant in both rodents and laurasiatherians:")
print(intersect(rod_sig_dec$X, lau_sig_dec$X))
print("Genes significant in both rodents and primates:")
print(intersect(rod_sig_dec$X, pri_sig_dec$X))
print("Genes significant in both primates and laurasiatherians:")
print(intersect(pri_sig_dec$X, lau_sig_dec$X))


#Get gene subsets
load("/ihome/nclark/emk270/MAMMAL_RTM_RER/GOdata_ensemblHs.rds")
GO_spermatogenesis<-unlist(genesInTerm(myGOdata, "GO:0007283"))
humanAtlasTestisSpecific<-read.table("/ihome/nclark/emk270/MAMMAL_RTM_RER/tissue_category_rna_testis_DetectedOnlyInTestis.tsv", sep="\t", header=TRUE)
GO_spermatogenesis_ts<-intersect(GO_spermatogenesis, humanAtlasTestisSpecific$Ensembl)

eut_sig_dec_GOsperm<-eut_sig_dec[which(eut_sig_dec$X %in% GO_spermatogenesis),]
pri_sig_dec_GOsperm<-pri_sig_dec[which(pri_sig_dec$X %in% GO_spermatogenesis),]
rod_sig_dec_GOsperm<-rod_sig_dec[which(rod_sig_dec$X %in% GO_spermatogenesis),]
lau_sig_dec_GOsperm<-lau_sig_dec[which(lau_sig_dec$X %in% GO_spermatogenesis),]

eut_sig_dec_GOsperm_ts<-eut_sig_dec[which(eut_sig_dec$X %in% GO_spermatogenesis_ts),]
pri_sig_dec_GOsperm_ts<-pri_sig_dec[which(pri_sig_dec$X %in% GO_spermatogenesis_ts),]
rod_sig_dec_GOsperm_ts<-rod_sig_dec[which(rod_sig_dec$X %in% GO_spermatogenesis_ts),]
lau_sig_dec_GOsperm_ts<-lau_sig_dec[which(lau_sig_dec$X %in% GO_spermatogenesis_ts),]

#Make Venn Diagrams for GO spermatogenesis and testis-specific datasets
catnames<-c("Eutherian", "Primate", "Rodent", "Laurasiatherian")
#vd<-venn.diagram(x=list(eut_sig_dec_GOsperm$X, pri_sig_dec_GOsperm$X, rod_sig_dec_GOsperm$X, lau_sig_dec_GOsperm$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.GOspermSigDec.all241datasets.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")
vd<-venn.diagram(x=list(eut_sig_dec_GOsperm$X, pri_sig_dec_GOsperm$X, rod_sig_dec_GOsperm$X, lau_sig_dec_GOsperm$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.GOspermSigDec.all241datasets.keepEarlyStops.noTraitWinsor.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")


catnames<-c("Primate", "Rodent", "Laurasiatherian")
#vd<-venn.diagram(x=list(eut_sig_dec_GOsperm$X, pri_sig_dec_GOsperm$X, rod_sig_dec_GOsperm$X, lau_sig_dec_GOsperm$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.GOspermSigDec.all241datasets.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")
vd<-venn.diagram(x=list(pri_sig_dec_GOsperm$X, rod_sig_dec_GOsperm$X, lau_sig_dec_GOsperm$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5), filename="./enrich_vis/venn_diagram.GOspermSigDec.all241datasets_NO_EUTH.keepEarlyStops.noTraitWinsor.png", imagetype="png", output=TRUE, fill=c("#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")




catnames<-c("Eutherian", "Primate", "Rodent", "Laurasiatherian")
#vd<-venn.diagram(x=list(eut_sig_dec_GOsperm_ts$X, pri_sig_dec_GOsperm_ts$X, rod_sig_dec_GOsperm_ts$X, lau_sig_dec_GOsperm_ts$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.GOspermTSSigDec.all241datasets.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")
vd<-venn.diagram(x=list(eut_sig_dec_GOsperm_ts$X, pri_sig_dec_GOsperm_ts$X, rod_sig_dec_GOsperm_ts$X, lau_sig_dec_GOsperm_ts$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.GOspermTSSigDec.all241datasets.keepEarlyStops.noTraitWinsor.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")

catnames<-c("Primate", "Rodent", "Laurasiatherian")
#vd<-venn.diagram(x=list(eut_sig_dec_GOsperm_ts$X, pri_sig_dec_GOsperm_ts$X, rod_sig_dec_GOsperm_ts$X, lau_sig_dec_GOsperm_ts$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.GOspermTSSigDec.all241datasets.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")
vd<-venn.diagram(x=list(pri_sig_dec_GOsperm_ts$X, rod_sig_dec_GOsperm_ts$X, lau_sig_dec_GOsperm_ts$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5), filename="./enrich_vis/venn_diagram.GOspermTSSigDec.all241datasets_NO_EUTH.keepEarlyStops.noTraitWinsor.png", imagetype="png", output=TRUE, fill=c("#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")



print("GO spermatogenesis genes significant in both eutherians and primates:")
print(intersect(eut_sig_dec_GOsperm$X, pri_sig_dec_GOsperm$X))
print("GO spermatogenesis genes significant in both eutherians and rodents:")
print(intersect(eut_sig_dec_GOsperm$X, rod_sig_dec_GOsperm$X))
print("GO spermatogenesis genes significant in both eutherians and laurasiatherians:")
print(intersect(eut_sig_dec_GOsperm$X, lau_sig_dec_GOsperm$X))
print("GO spermatogenesis genes significant in both rodents and laurasiatherians:")
print(intersect(rod_sig_dec_GOsperm$X, lau_sig_dec_GOsperm$X))
print("GO spermatogenesis genes significant in both rodents and primates")
print(intersect(rod_sig_dec_GOsperm$X, pri_sig_dec_GOsperm$X))
print("GO spermatogenesis genes significant in both primates and laurasiatherians:")
print(intersect(pri_sig_dec_GOsperm$X, lau_sig_dec_GOsperm$X))

#Test for correlations of Rho values among spermatogenesis genes
eut_GOsperm<-eut[which(eut$X %in% GO_spermatogenesis),]
pri_GOsperm<-pri[which(pri$X %in% GO_spermatogenesis),]
rod_GOsperm<-rod[which(rod$X %in% GO_spermatogenesis),]
lau_GOsperm<-lau[which(lau$X %in% GO_spermatogenesis),]

keep<-Reduce(intersect, list(eut_GOsperm$X, pri_GOsperm$X, rod_GOsperm$X, lau_GOsperm$X))

eut_GOsperm_keep<-eut_GOsperm[which(eut_GOsperm$X %in% keep),]
eut_GOsperm_order<-eut_GOsperm_keep[order(eut_GOsperm_keep$X),]
pri_GOsperm_keep<-pri_GOsperm[which(pri_GOsperm$X %in% keep),]
pri_GOsperm_order<-pri_GOsperm_keep[order(pri_GOsperm_keep$X),]
rod_GOsperm_keep<-rod_GOsperm[which(rod_GOsperm$X %in% keep),]
rod_GOsperm_order<-rod_GOsperm_keep[order(rod_GOsperm_keep$X),]
lau_GOsperm_keep<-lau_GOsperm[which(lau_GOsperm$X %in% keep),]
lau_GOsperm_order<-lau_GOsperm_keep[order(lau_GOsperm_keep$X),]

stopifnot(all.equal(eut_GOsperm_order$X, pri_GOsperm_order$X))
stopifnot(all.equal(eut_GOsperm_order$X, rod_GOsperm_order$X))
stopifnot(all.equal(eut_GOsperm_order$X, lau_GOsperm_order$X))

eutVSpri_GOsperm<-cor.test(x=eut_GOsperm_order$Rho, y=pri_GOsperm_order$Rho)
eutVSrod_GOsperm<-cor.test(x=eut_GOsperm_order$Rho, y=rod_GOsperm_order$Rho)
eutVSlau_GOsperm<-cor.test(x=eut_GOsperm_order$Rho, y=lau_GOsperm_order$Rho)
priVSrod_GOsperm<-cor.test(x=pri_GOsperm_order$Rho, rod_GOsperm_order$Rho)
priVSlau_GOsperm<-cor.test(x=pri_GOsperm_order$Rho, lau_GOsperm_order$Rho)
rodVSlau_GOsperm<-cor.test(x=rod_GOsperm_order$Rho, lau_GOsperm_order$Rho)

print(eutVSpri_GOsperm)
print(eutVSrod_GOsperm)
print(eutVSlau_GOsperm)
print(priVSrod_GOsperm)
print(priVSlau_GOsperm)
print(rodVSlau_GOsperm)

#pdf("correlations.GOsperm.all241datasets.pdf", onefile=TRUE)
pdf("correlations.GOsperm.all241datasets.keepEarlyStops.noTraitWinsor.pdf", onefile=TRUE)
plot(x=eut_GOsperm_order$Rho, y=pri_GOsperm_order$Rho, main="All eutherians vs primates\nGO spermatogenesis genes")
text(x=-0.25, y=0.4, labels=paste("Cor =", signif(eutVSpri_GOsperm$estimate, 3)))
text(x=-0.25, y=0.3, labels=paste("P =", signif(eutVSpri_GOsperm$p.value, 3)))
plot(x=eut_GOsperm_order$Rho, y=rod_GOsperm_order$Rho, main="All eutherians vs rodents\nGO spermatogenesis genes")
text(x=-0.25, y=0.6, labels=paste("Cor =", signif(eutVSrod_GOsperm$estimate,3)))
text(x=-0.25, y=0.5, labels=paste("P =", signif(eutVSrod_GOsperm$p.value,3)))
plot(x=eut_GOsperm_order$Rho, y=lau_GOsperm_order$Rho, main="All eutherians vs laurasiatherians\nGO spermatogenesis genes")
text(x=-0.25, y=0.4, labels=paste("Cor =", signif(eutVSlau_GOsperm$estimate,3)))
text(x=-0.25, y=0.3, labels=paste("P =", signif(eutVSlau_GOsperm$p.value,3)))
plot(x=pri_GOsperm_order$Rho, y=rod_GOsperm_order$Rho, main="Primates vs rodents\nGO spermatogenesis genes")
text(x=-0.5, y=0.6, labels=paste("Cor =", signif(priVSrod_GOsperm$estimate,3)))
text(x=-0.5, y=0.5, labels=paste("P =", signif(priVSrod_GOsperm$p.value,3)))
plot(x=pri_GOsperm_order$Rho, y=lau_GOsperm_order$Rho, main="Primates vs laurasiatherians\nGO spermatogenesis genes")
text(x=-0.5, y=0.4, labels=paste("Cor =", signif(priVSlau_GOsperm$estimate,3)))
text(x=-0.5, y=0.3, labels=paste("P =", signif(priVSlau_GOsperm$p.value,3)))
plot(x=rod_GOsperm_order$Rho, y=lau_GOsperm_order$Rho, main="Rodents vs laurasiatherians\nGO spermatogenesis genes")
text(x=-0.5, y=0.4, labels=paste("Cor =", signif(rodVSlau_GOsperm$estimate,3)))
text(x=-0.5, y=0.3, labels=paste("P =", signif(rodVSlau_GOsperm$p.value,3)))
dev.off()

#Compare stat (signed log10 P-vals) across taxa
print("Testing for correlations among sign log10 P-values")
rhosign<-sign(eut_GOsperm_order$Rho)
eut_GOsperm_order$signLogP<-(-rhosign*log10(eut_GOsperm_order$P))
rhosign<-sign(pri_GOsperm_order$Rho)
pri_GOsperm_order$signLogP<-(-rhosign*log10(pri_GOsperm_order$P))
rhosign<-sign(rod_GOsperm_order$Rho)
rod_GOsperm_order$signLogP<-(-rhosign*log10(rod_GOsperm_order$P))
rhosign<-sign(lau_GOsperm_order$Rho)
lau_GOsperm_order$signLogP<-(-rhosign*log10(lau_GOsperm_order$P))

eutVSpri_GOsperm_stat<-cor.test(x=eut_GOsperm_order$signLogP, y=pri_GOsperm_order$signLogP)
eutVSrod_GOsperm_stat<-cor.test(x=eut_GOsperm_order$signLogP, y=rod_GOsperm_order$signLogP)
eutVSlau_GOsperm_stat<-cor.test(x=eut_GOsperm_order$signLogP, y=lau_GOsperm_order$signLogP)
priVSrod_GOsperm_stat<-cor.test(x=pri_GOsperm_order$signLogP, rod_GOsperm_order$signLogP)
priVSlau_GOsperm_stat<-cor.test(x=pri_GOsperm_order$signLogP, lau_GOsperm_order$signLogP)
rodVSlau_GOsperm_stat<-cor.test(x=rod_GOsperm_order$signLogP, lau_GOsperm_order$signLogP)

print(eutVSpri_GOsperm_stat)
print(eutVSrod_GOsperm_stat)
print(eutVSlau_GOsperm_stat)
print(priVSrod_GOsperm_stat)
print(priVSlau_GOsperm_stat)
print(rodVSlau_GOsperm_stat)

#pdf("correlations.GOsperm.all241datasets.pdf", onefile=TRUE)
pdf("correlations.GOsperm.signLogP.all241datasets.keepEarlyStops.noTraitWinsor.pdf", onefile=TRUE)
plot(x=eut_GOsperm_order$signLogP, y=pri_GOsperm_order$signLogP, main="All eutherians vs primates\nGO spermatogenesis genes")
text(x=-2.5, y=1.1, labels=paste("Cor =", signif(eutVSpri_GOsperm_stat$estimate, 3)))
text(x=-2.5, y=0.9, labels=paste("P =", signif(eutVSpri_GOsperm_stat$p.value, 3)))
plot(x=eut_GOsperm_order$signLogP, y=rod_GOsperm_order$signLogP, main="All eutherians vs rodents\nGO spermatogenesis genes")
text(x=-2.5, y=3, labels=paste("Cor =", signif(eutVSrod_GOsperm_stat$estimate,3)))
text(x=-2.5, y=2.5, labels=paste("P =", signif(eutVSrod_GOsperm_stat$p.value,3)))
plot(x=eut_GOsperm_order$signLogP, y=lau_GOsperm_order$signLogP, main="All eutherians vs laurasiatherians\nGO spermatogenesis genes")
text(x=-2.5, y=2, labels=paste("Cor =", signif(eutVSlau_GOsperm_stat$estimate,3)))
text(x=-2.5, y=1.5, labels=paste("P =", signif(eutVSlau_GOsperm_stat$p.value,3)))
plot(x=pri_GOsperm_order$signLogP, y=rod_GOsperm_order$signLogP, main="Primates vs rodents\nGO spermatogenesis genes")
text(x=-2.5, y=3, labels=paste("Cor =", signif(priVSrod_GOsperm_stat$estimate,3)))
text(x=-2.5, y=2.5, labels=paste("P =", signif(priVSrod_GOsperm_stat$p.value,3)))
plot(x=pri_GOsperm_order$signLogP, y=lau_GOsperm_order$signLogP, main="Primates vs laurasiatherians\nGO spermatogenesis genes")
text(x=-2.5, y=2.1, labels=paste("Cor =", signif(priVSlau_GOsperm_stat$estimate,3)))
text(x=-2.5, y=1.9, labels=paste("P =", signif(priVSlau_GOsperm_stat$p.value,3)))
plot(x=rod_GOsperm_order$signLogP, y=lau_GOsperm_order$signLogP, main="Rodents vs laurasiatherians\nGO spermatogenesis genes")
text(x=-2.5, y=2.1, labels=paste("Cor =", signif(rodVSlau_GOsperm_stat$estimate,3)))
text(x=-2.5, y=1.9, labels=paste("P =", signif(rodVSlau_GOsperm_stat$p.value,3)))
dev.off()

#Compare 250 most decelerated genes across all datasets
eut_sig_dec_sort<-eut_sig_dec[order(eut_sig_dec$P),]
pri_sig_dec_sort<-pri_sig_dec[order(pri_sig_dec$P),]
rod_sig_dec_sort<-rod_sig_dec[order(rod_sig_dec$P),]
lau_sig_dec_sort<-lau_sig_dec[order(lau_sig_dec$P),]

eut_sig_dec_250<-eut_sig_dec_sort[1:250,]
pri_sig_dec_250<-pri_sig_dec_sort[1:250,]
rod_sig_dec_250<-rod_sig_dec_sort[1:250,]
lau_sig_dec_250<-lau_sig_dec_sort[1:250,]

catnames<-c("Eutherian", "Primate", "Rodent", "Laurasiatherian")
#vd<-venn.diagram(x=list(eut_sig_dec_250$X, pri_sig_dec_250$X, rod_sig_dec_250$X, lau_sig_dec_250$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.top250SigDec.all241datasets.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")
vd<-venn.diagram(x=list(eut_sig_dec_250$X, pri_sig_dec_250$X, rod_sig_dec_250$X, lau_sig_dec_250$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.top250SigDec.all241datasets.keepEarlyStops.noTraitWinsor.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")



#compare top 50 and 100
#Compare 250 most decelerated genes across all datasets
eut_sig_dec_sort<-eut_sig_dec[order(eut_sig_dec$P),]
pri_sig_dec_sort<-pri_sig_dec[order(pri_sig_dec$P),]
rod_sig_dec_sort<-rod_sig_dec[order(rod_sig_dec$P),]
lau_sig_dec_sort<-lau_sig_dec[order(lau_sig_dec$P),]

eut_sig_dec_100<-eut_sig_dec_sort[1:100,]
pri_sig_dec_100<-pri_sig_dec_sort[1:100,]
rod_sig_dec_100<-rod_sig_dec_sort[1:100,]
lau_sig_dec_100<-lau_sig_dec_sort[1:100,]

catnames<-c("Eutherian", "Primate", "Rodent", "Laurasiatherian")
#vd<-venn.diagram(x=list(eut_sig_dec_250$X, pri_sig_dec_250$X, rod_sig_dec_250$X, lau_sig_dec_250$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.top250SigDec.all241datasets.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")
vd<-venn.diagram(x=list(eut_sig_dec_100$X, pri_sig_dec_100$X, rod_sig_dec_100$X, lau_sig_dec_100$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="./enrich_vis/venn_diagram.top100SigDec.all241datasets.keepEarlyStops.noTraitWinsor.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")


#Compare 250 most decelerated genes across all datasets
eut_sig_dec_sort<-eut_sig_dec[order(eut_sig_dec$P),]
pri_sig_dec_sort<-pri_sig_dec[order(pri_sig_dec$P),]
rod_sig_dec_sort<-rod_sig_dec[order(rod_sig_dec$P),]
lau_sig_dec_sort<-lau_sig_dec[order(lau_sig_dec$P),]

eut_sig_dec_50<-eut_sig_dec_sort[1:50,]
pri_sig_dec_50<-pri_sig_dec_sort[1:50,]
rod_sig_dec_50<-rod_sig_dec_sort[1:50,]
lau_sig_dec_50<-lau_sig_dec_sort[1:50,]

catnames<-c("Eutherian", "Primate", "Rodent", "Laurasiatherian")
#vd<-venn.diagram(x=list(eut_sig_dec_250$X, pri_sig_dec_250$X, rod_sig_dec_250$X, lau_sig_dec_250$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="venn_diagram.top250SigDec.all241datasets.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")
vd<-venn.diagram(x=list(eut_sig_dec_50$X, pri_sig_dec_50$X, rod_sig_dec_50$X, lau_sig_dec_50$X), category.names=catnames, cat.dist=c(0.5,0.5,0.5,0.5), filename="./enrich_vis/venn_diagram.top50SigDec.all241datasets.keepEarlyStops.noTraitWinsor.png", imagetype="png", output=TRUE, fill=c("#332288", "#DDCC77", "#CC6677", "#44AA99"), alpha=c(0.3, 0.3, 0.3, 0.3), cex=0.9, fontfamily="sans")











#Make a table of genes significant in at least one clade - which clades are they sig in? Are they GO spermatogenesis genes?
sigGenes<-Reduce(union, list(eut_sig_dec$X, pri_sig_dec$X, rod_sig_dec$X, lau_sig_dec$X))
sigEut<-c()
sigPri<-c()
sigRod<-c()
sigLau<-c()
spermGO<-c()
for(i in sigGenes){
        if(i %in% eut_sig_dec$X){
                sigEut<-c(sigEut, "yes")
        } else{
                sigEut<-c(sigEut, "no")
        }
        if(i %in% pri_sig_dec$X){
                sigPri<-c(sigPri, "yes")
        } else{
                sigPri<-c(sigPri, "no")
        }
        if(i %in% rod_sig_dec$X){
                sigRod<-c(sigRod, "yes")
        } else{
                sigRod<-c(sigRod, "no")
        }
        if(i %in% lau_sig_dec$X){
                sigLau<-c(sigLau, "yes")
        } else{
                sigLau<-c(sigLau, "no")
        }
	if(i %in% GO_spermatogenesis){
		spermGO<-c(spermGO, "yes")
	} else{
		spermGO<-c(spermGO, "no")
	}
}
sigGenesTable<-as.data.frame(cbind(sigGenes, sigEut, sigPri, sigRod, sigLau, spermGO))
write.csv(sigGenesTable, file="genesSigInAtLeastOneClade.keepEarlyStops.noTraitWinsor.csv", quote=FALSE, row.names=FALSE)

#GO enrichment for genes sig in both Eutherians and Laurasiatherians, compared to background set of genes in both dataset
bg<-intersect(eut$X, lau$X)
fg<-intersect(eut_sig_dec$X, lau_sig_dec$X)
inSubset<-c()
for(gene in bg){
        if(gene %in% fg){
                inSubset<-c(inSubset, 1)
        } else{
                inSubset<-c(inSubset, 0)
        }
}
if(length(fg) != length(which(inSubset==1))){
        print("ERROR: length of 'subset_dec' file and 'inSubset' vector don't match up!")
}
names(inSubset)<-bg
geneSelFunc<-function(iO){
        return(iO==1)
}
myGOdata<-new("topGOdata", description="Plt05Vtotal_dec", ontology="BP", allGenes=inSubset, geneSel=geneSelFunc, annot=annFUN.org, mapping="org.Hs.eg.db", ID="ensembl")
result<-runTest(myGOdata, algorithm = "classic", statistic = "fisher")
resultData<-GenTable(myGOdata, classFisher=result, orderBy="classFisher", ranksOf="classicFisher", topNodes=length(score(result)))
resultData$p.adj<-p.adjust(resultData$classFisher, method="BH")
print(head(resultData))
write.csv(resultData, file=paste0("RERconverge_GOenrichment.sigInEuthAndLaura.translated.keepEarlyStops.noTraitWinsor.csv"), row.names=FALSE)

