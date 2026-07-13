#PURPOSE: Generate plots of enrichment density by RERconverge rank -- can plot with either RERconverge Rho or Stat statistic

library(topGO)
library(ggplot2)
source("/ix3/nclark/ekopania/MURINAE_REVISIONS/RERconverge_noParalogs/Mariabarcodeplotfinal.R")

#rhoMat<-c()
signedPMat<-c()

#res<-read.csv("241continuousGENES_RERconvergeResults.csv", header=TRUE)
#res0<-read.csv("translated241_results_log_rts_csv_trusted_241continuous.csv", header=TRUE)
#res0<-read.csv("241continuous_translated_RER_results.keepEarlyStops.csv", header=TRUE)
#res0<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/RTS_RER_translated_keepEarlyStops_noAsk_results.csv", header=TRUE)
res0 = read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/binary.genetic.mating.systems_nocomp.as.foreground_results.csv", header=TRUE)
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
#annots<-annFUN.org("BP", mapping="org.Hs.eg.db", ID="ensembl")
#Load GO data from previous topGO runs
load("/ihome/nclark/emk270/MAMMAL_RTM_RER/GOdata_ensemblHs.rds")
GO_spermatogenesis<-unlist(genesInTerm(myGOdata, "GO:0007283"))
GO_gameteGeneration<-unlist(genesInTerm(myGOdata, "GO:0007276"))
GO_acrosomeReaction<-unlist(genesInTerm(myGOdata, "GO:0007340"))
GO_prematureAcrosomeLoss<-unlist(genesInTerm(myGOdata, "GO:0061948"))
GO_fertilization<-unlist(genesInTerm(myGOdata, "GO:0009566"))
GO_ovulationCycleProcess<-unlist(genesInTerm(myGOdata, "GO:0022602"))
GO_spermMotility<-unlist(genesInTerm(myGOdata, "GO:0097722"))
humanAtlasTestisSpecific<-read.table("/ihome/nclark/emk270/MAMMAL_RTM_RER/tissue_category_rna_testis_DetectedOnlyInTestis.tsv", sep="\t", header=TRUE)
GO_spermatogenesis_ts<-intersect(GO_spermatogenesis, humanAtlasTestisSpecific$Ensembl)
#Others to try:
#ion channels? Ask Justine about other hits

print("Running plot function...")
pdf("./enrich_vis/enrichment_density.241.BINARY.nocompforeground.translated.KeepEarlyStop.noAsk.span0_2.byStat.noTraitWinsor.pdf", height=5, width=7, onefile=TRUE)
makeMultiGSEAplot(signedPMat, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis")
makeMultiGSEAplot(signedPMat, GO_gameteGeneration, span=.2, title="GO category: gamete generation")
makeMultiGSEAplot(signedPMat, GO_acrosomeReaction, span=.2, title="GO category: acrosome reaction")
makeMultiGSEAplot(signedPMat, GO_prematureAcrosomeLoss, span=.2, title="GO category: premature acrosome loss")
makeMultiGSEAplot(signedPMat, GO_fertilization, span=.2, title="GO category: fertilization")
makeMultiGSEAplot(signedPMat, GO_ovulationCycleProcess, span=.2, title="GO category: ovulation cycle process")
makeMultiGSEAplot(signedPMat, GO_spermMotility, span=.2, title="GO category: sperm motility")
makeMultiGSEAplot(signedPMat, GO_spermatogenesis_ts, span=.2, title="GO category: spermatogenesis AND testis-specific")
dev.off()

#Repeat for primates
print("Repeat for primates (Ports paper)...")
res_ports0<-read.csv("/ihome/nclark/emk270/MAMMAL_RTM_RER/portsAndJensenSeaman_RERconvergeResults_suppTable3.withGeneID.csv", header=TRUE)
#Remove NAs
res_ports<-res_ports0[which(!(is.na(res_ports0$Rho))),]
rhosign<-sign(res_ports$Rho)
res_ports$signLogP<-(-rhosign*log10(res_ports$P))

signedPMat_ports<-c()
signedPMat_ports<-cbind(signedPMat_ports, res_ports$signLogP)
rownames(signedPMat_ports)<-res_ports$geneID

print("signedPMat_ports info:")
colnames(signedPMat_ports)<-c("Ports")
print(paste("Matrix?", is.matrix(signedPMat_ports)))
print(dim(signedPMat_ports))
print(head(signedPMat_ports))
pdf("./enrich_vis/enrichment_density.241continuous.Ports_primates.span0_2.byStat.pdf", height=5, width=7, onefile=TRUE)
makeMultiGSEAplot(signedPMat_ports, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis; Ports primate data")
makeMultiGSEAplot(signedPMat_ports, GO_spermatogenesis_ts, span=.2, title="GO category: spermatogenesis AND testis-specific\nPorts primate data")
dev.off()

#Repeat for murines
print("Repeat for murines (Kopania et al.)...")
res_kopania0<-read.csv("/ihome/nclark/emk270/MAMMAL_RTM_RER/kopaniaEtAl_ST9_molec_evo_all.withGeneID.csv", header=TRUE)
#Remove NAs
res_kopania<-res_kopania0[which(!(is.na(res_kopania0$RERconverge_Rho))),]
rhosign<-sign(res_kopania$RERconverge_Rho)
res_kopania$RERconverge_signLogP<-(-rhosign*log10(res_kopania$RERconverge_P))

signedPMat_kopania<-c()
#Negative because I set small RTM species as foreground
signedPMat_kopania<-cbind(signedPMat_kopania, -res_kopania$RERconverge_signLogP)
rownames(signedPMat_kopania)<-res_kopania$humID

print("signedPMat_kopania info:")
colnames(signedPMat_kopania)<-c("Kopania")
print(paste("Matrix?", is.matrix(signedPMat_kopania)))
print(dim(signedPMat_kopania))
print(head(signedPMat_kopania))
pdf("./enrich_vis/enrichment_density.241continuous.Kopania_murines.span0_2.byStat.pdf", height=5, width=7, onefile=TRUE)
makeMultiGSEAplot(signedPMat_kopania, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis; Kopania murine data")
makeMultiGSEAplot(signedPMat_kopania, GO_spermatogenesis_ts, span=.2, title="GO category: spermatogenesis AND testis-specific\nKopania murine data")
dev.off()

#Combine all 3 in one plot
print("Combine all 3 datasets in one plot...")

keep<-Reduce(intersect, list(res$X, res_ports$geneID, res_kopania$humID))
print(paste("Number of genes in all 3 datasets:", length(keep)))

res_filtered<-res[which(res$X %in% keep),]
signedPMat<-c()
signedPMat<-cbind(signedPMat, res_filtered$signLogP)
rownames(signedPMat)<-res_filtered$X
#print("signLogPMat:")
#print(head(signedPMat))

res_ports_filtered<-res_ports[which(res_ports$geneID %in% keep),]
res_ports_sorted<-res_ports_filtered[match(rownames(signedPMat), res_ports_filtered$geneID),]
#print("res_ports_sorted:")
#print(head(res_ports_sorted))
stopifnot(all.equal(rownames(signedPMat), res_ports_sorted$geneID))
signedPMat<-cbind(signedPMat,res_ports_sorted$signLogP)

res_kopania_filtered<-res_kopania[which(res_kopania$humID %in% keep),]
res_kopania_sorted<-res_kopania_filtered[match(rownames(signedPMat), res_kopania_filtered$humID),]
#print("res_kopania_sorted")
#print(head(res_kopania_sorted))
stopifnot(all.equal(rownames(signedPMat), res_kopania_sorted$humID))
#Negative because I set small RTS species as foreground
signedPMat<-cbind(signedPMat, -res_kopania_sorted$RERconverge_signLogP)

print("signedPMat info:")
colnames(signedPMat)<-c("eutherians", "Ports_primates", "Kopania_murines")
print(paste("Matrix?", is.matrix(signedPMat)))
print(dim(signedPMat))
print(head(signedPMat))
pdf("./enrich_vis/enrichment_density.241continuous.datasets_combined.span0_2.byStat.keepEarlyStops.noTraitWinsor.pdf", height=5, width=7, onefile=TRUE)
makeMultiGSEAplot(signedPMat, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis", cols=c("#332288","#88CCEE","#CC6677"))
makeMultiGSEAplot(signedPMat, GO_spermatogenesis_ts, span=.2, title="GO category: spermatogenesis AND testis-specific",cols=c("#332288","#88CCEE","#CC6677"))
#Plot separately but with only the genes present in all 3 datasets
signedPMat_euth<-as.matrix(signedPMat[,"eutherians"])
colnames(signedPMat_euth)<-c("eutherians")
print(head(signedPMat_euth))
makeMultiGSEAplot(signedPMat_euth, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis; eutherians only", cols=c("#332288"))
makeMultiGSEAplot(signedPMat_euth, GO_spermatogenesis_ts, span=.2, title="GO category: spermatogenesis AND testis-specific\neutherians only", cols=c("#332288"))
signedPMat_pri<-as.matrix(signedPMat[,"Ports_primates"])
colnames(signedPMat_pri)<-c("Ports_primates")
makeMultiGSEAplot(signedPMat_pri, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis; primates only",cols=c("#DDCC77"))
makeMultiGSEAplot(signedPMat_pri, GO_spermatogenesis_ts, span=.2, title="GO category: spermatogenesis AND testis-specific\nprimates only", cols=c("#DDCC77"))
signedPMat_mur<-as.matrix(signedPMat[,"Kopania_murines"])
colnames(signedPMat_mur)<-c("Kopania_murines")
makeMultiGSEAplot(signedPMat_mur, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis; murines only", cols=c("#CC6677"))
makeMultiGSEAplot(signedPMat_mur, GO_spermatogenesis_ts, span=.2, title="GO category: spermatogenesis AND testis-specific\nmurines only", cols=c("#CC6677"))
dev.off()

#Combine all eutherians with clade-specific results from zoonomia 241 dataset
#zoonomia_pri_res0<-read.csv("241CONTINUOUS_PRIMATES_ONLY/translated_43species_breadTree_LOGGED_MAREL_primate.csv", header=TRUE)
#zoonomia_pri_res<-zoonomia_pri_res0[which(!(is.na(zoonomia_pri_res0$stat))),]
#zoonomia_pri_res$signLogP<-zoonomia_pri_res$stat
zoonomia_pri_res0<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/translated_breadTree_primates_only_KeepEarlyStops_noAsk_RERresults.csv", header=TRUE)
zoonomia_pri_res<-zoonomia_pri_res0[which(!(is.na(zoonomia_pri_res0$Rho))),]
rhosign<-sign(zoonomia_pri_res$Rho)
zoonomia_pri_res$signLogP<-(-rhosign*log10(zoonomia_pri_res$P))

#zoonomia_rod_res0<-read.csv("241CONTINUOUS_RODENTS_ONLY/translated_breadTree_53RODENT_RERresults.csv", header=TRUE)
#zoonomia_rod_res<-zoonomia_rod_res0[which(!(is.na(zoonomia_rod_res0$stat))),]
#zoonomia_rod_res$signLogP<-zoonomia_rod_res$stat
zoonomia_rod_res0<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/translated_breadTree_rodents_only_KeepEarlyStops_noAsk_RERresults.csv", header=TRUE)
zoonomia_rod_res<-zoonomia_rod_res0[which(!(is.na(zoonomia_rod_res0$Rho))),]
rhosign<-sign(zoonomia_rod_res$Rho)
zoonomia_rod_res$signLogP<-(-rhosign*log10(zoonomia_rod_res$P))

#zoonomia_lau_res0<-read.csv("241CONTINUOUS_LAURASIATHERIA_ONLY/translated_breadTree_lauras_only_RERresults.csv", header=TRUE)
#zoonomia_lau_res<-zoonomia_lau_res0[which(!(is.na(zoonomia_lau_res0$stat))),]
#zoonomia_lau_res$signLogP<-zoonomia_lau_res$stat
zoonomia_lau_res0<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/translated_breadTree_lauras_only_KeepEarlyStops_noAsk_RERresults.csv", header=TRUE)
zoonomia_lau_res<-zoonomia_lau_res0[which(!(is.na(zoonomia_lau_res0$Rho))),]
rhosign<-sign(zoonomia_lau_res$Rho)
zoonomia_lau_res$signLogP<-(-rhosign*log10(zoonomia_lau_res$P))

keep<-Reduce(intersect, list(res$X, zoonomia_pri_res$X, zoonomia_rod_res$X, zoonomia_lau_res$X))
print(paste("Number of genes in all 4 datasets:", length(keep)))

res_filtered<-res[which(res$X %in% keep),]
zoonomia_pri_res_filtered<-zoonomia_pri_res[which(zoonomia_pri_res$X %in% keep),]
zoonomia_rod_res_filtered<-zoonomia_rod_res[which(zoonomia_rod_res$X %in% keep),]
zoonomia_lau_res_filtered<-zoonomia_lau_res[which(zoonomia_lau_res$X %in% keep),]

res_sorted<-res_filtered[order(res_filtered$X),]
zoonomia_pri_res_sorted<-zoonomia_pri_res_filtered[order(zoonomia_pri_res_filtered$X),]
zoonomia_rod_res_sorted<-zoonomia_rod_res_filtered[order(zoonomia_rod_res_filtered$X),]
zoonomia_lau_res_sorted<-zoonomia_lau_res_filtered[order(zoonomia_lau_res_filtered$X),]

stopifnot(all.equal(res_sorted$X, zoonomia_pri_res_sorted$X))
stopifnot(all.equal(res_sorted$X, zoonomia_rod_res_sorted$X))
stopifnot(all.equal(res_sorted$X, zoonomia_lau_res_sorted$X))
signedPMat_zoonomia<-cbind(res_sorted$signLogP, zoonomia_pri_res_sorted$signLogP, zoonomia_rod_res_sorted$signLogP, zoonomia_lau_res_sorted$signLogP)
rownames(signedPMat_zoonomia)<-res_sorted$X

print("signedPMat info:")
colnames(signedPMat_zoonomia)<-c("eutherians241", "primates241", "rodents241", "lauras241")
print(paste("Matrix?", is.matrix(signedPMat_zoonomia)))
print(dim(signedPMat_zoonomia))
print(head(signedPMat_zoonomia))
pdf("./enrich_vis/enrichment_density.241continuous.separate_clades.span0_2.byStat.keepEarlyStops.noTraitWinsor.pdf", height=5, width=7, onefile=TRUE)
makeMultiGSEAplot(signedPMat_zoonomia, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis\nseparate clades w/ Zoonomia alignments", cols=c("#332288","#DDCC77","#882255","#44AA99"))
makeMultiGSEAplot(signedPMat_zoonomia, GO_spermatogenesis_ts, span=.2, title="GO category: spermatogenesis AND testis-specific\nseparate clades w/ Zoonomia alignments", cols=c("#332288","#DDCC77","#882255","#44AA99"))
dev.off()

#Combine all datasets (zoonomia 241 dataset all eutherians, primates, rodents, Laurasiatheria, Ports primats, Kopania murines)
keep<-Reduce(intersect, list(res$X, zoonomia_pri_res$X, zoonomia_rod_res$X, zoonomia_lau_res$X, res_ports$geneID, res_kopania$humID))
print(paste("Number of genes in all datasets:", length(keep)))

res_filtered<-res[which(res$X %in% keep),]
zoonomia_pri_res_filtered<-zoonomia_pri_res[which(zoonomia_pri_res$X %in% keep),]
zoonomia_rod_res_filtered<-zoonomia_rod_res[which(zoonomia_rod_res$X %in% keep),]
zoonomia_lau_res_filtered<-zoonomia_lau_res[which(zoonomia_lau_res$X %in% keep),]
ports_res_filtered<-res_ports[which(res_ports$geneID %in% keep),]
kopania_res_filtered<-res_kopania[which(res_kopania$humID %in% keep),]

res_sorted<-res_filtered[order(res_filtered$X),]
zoonomia_pri_res_sorted<-zoonomia_pri_res_filtered[match(res_sorted$X, zoonomia_pri_res_filtered$X),]
zoonomia_rod_res_sorted<-zoonomia_rod_res_filtered[match(res_sorted$X, zoonomia_rod_res_filtered$X),]
zoonomia_lau_res_sorted<-zoonomia_lau_res_filtered[match(res_sorted$X, zoonomia_lau_res_filtered$X),]
ports_res_sorted<-ports_res_filtered[match(res_sorted$X, ports_res_filtered$geneID),]
kopania_res_sorted<-kopania_res_filtered[match(res_sorted$X, kopania_res_filtered$humID),]

stopifnot(all.equal(res_sorted$X, zoonomia_pri_res_sorted$X))
stopifnot(all.equal(res_sorted$X, zoonomia_rod_res_sorted$X))
stopifnot(all.equal(res_sorted$X, zoonomia_lau_res_sorted$X))
stopifnot(all.equal(res_sorted$X, ports_res_sorted$geneID))
stopifnot(all.equal(res_sorted$X, kopania_res_sorted$humID))
signedPMat_all<-cbind(res_sorted$signLogP, zoonomia_pri_res_sorted$signLogP, zoonomia_rod_res_sorted$signLogP, zoonomia_lau_res_sorted$signLogP, ports_res_sorted$signLogP, -kopania_res_sorted$RERconverge_signLogP)
rownames(signedPMat_all)<-res_sorted$X

print("signedPMat info:")
colnames(signedPMat_all)<-c("eutherians241", "primates241", "rodents241", "lauras241", "primates_Ports", "murines_Kopania")
print(paste("Matrix?", is.matrix(signedPMat_all)))
print(dim(signedPMat_all))
print(head(signedPMat_all))
pdf("./enrich_vis/enrichment_density.241continuous.all_datasets.span0_2.byStat.keepEarlyStops.noTraitWinsor.pdf", height=5, width=7, onefile=TRUE)
makeMultiGSEAplot(signedPMat_all, GO_spermatogenesis, span=.2, title="GO category: spermatogenesis\nall datasets", cols=c("#332288","#DDCC77","#882255","#44AA99", "#88CCEE","#CC6677"))
makeMultiGSEAplot(signedPMat_all, GO_spermatogenesis_ts, span=.2, title="GO category: spermatogenesis AND testis-specific\nall datasets", cols=c("#332288","#DDCC77","#882255","#44AA99", "#88CCEE","#CC6677"))
dev.off()
