#PURPOSE: run topGO enrichment functions on RERconverge output files to determine enrichment across runs
library(topGO)
#res<-read.csv("241continuous_translated_RER_results.keepEarlyStops.csv", header=TRUE)
#res = read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/binary.genetic.mating.systems_nocomp.as.foreground_results.csv")
#res = read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/binary.genetic.mating.systems_COMP.as.foreground_results.csv")
#res = read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/binary.genetic.mating.systems.LESS.STRICT_nonCOMP.as.foreground_results.csv")

#res<-read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/RTS_RER_translated_keepEarlyStops_noAsk_results.csv", header=TRUE)
#res = read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/translated_breadTree_lauras_only_KeepEarlyStops_noAsk_RERresults.csv", header=TRUE)
#res = read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/translated_breadTree_primates_only_KeepEarlyStops_noAsk_RERresults.csv",header=TRUE)
res = read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/translated_breadTree_rodents_only_KeepEarlyStops_noAsk_RERresults.csv",header=TRUE)



print(head(res))

#Get top 250 most accelerated, based on significance
res_acc<-res[which(res$Rho > 0),]
print(dim(res_acc))
#order_acc<-res_acc[order(res_acc$p.adj),]
order_acc<-res_acc[order(res_acc$P),]
print(head(order_acc))
print(dim(order_acc))
print(order_acc[240:250,])
subset_acc<-order_acc$X[1:250]
#subset_acc<-order_acc$X[1:100]
#subset_acc<-order_acc$X[1:50]
print(head(subset_acc))
print(tail(subset_acc))

#Get top 250 most decelerated, based on significance
res_dec<-res[which(res$Rho < 0),]
print(dim(res_dec))
#order_dec<-res_dec[order(res_dec$p.adj),]
order_dec<-res_dec[order(res_dec$P),]
print(head(order_dec))
print(dim(order_dec))
print(order_dec[240:250,])
subset_dec<-order_dec$X[1:250]
#subset_dec<-order_dec$X[1:100]
#subset_dec<-order_dec$X[1:50]
print(head(subset_dec))
print(tail(subset_dec))

#Accelerated only
total<-order_acc$X
inSubset<-c()
for(gene in total){
        if(gene %in% subset_acc){
                inSubset<-c(inSubset, 1)
        } else{
                inSubset<-c(inSubset, 0)
        }
}
if(length(subset_acc) != length(which(inSubset==1))){
        print("ERROR: length of 'subset_acc' file and 'inSubset' vector don't match up!")
}
names(inSubset)<-total
geneSelFunc<-function(iO){
        return(iO==1)
}
myGOdata<-new("topGOdata", description="Plt05Vtotal_acc", ontology="BP", allGenes=inSubset, geneSel=geneSelFunc, annot=annFUN.org, mapping="org.Hs.eg.db", ID="ensembl")
result<-runTest(myGOdata, algorithm = "classic", statistic = "fisher")
resultData<-GenTable(myGOdata, classFisher=result, orderBy="classFisher", ranksOf="classicFisher", topNodes=length(score(result)))
resultData$p.adj<-p.adjust(resultData$classFisher, method="BH")
print(head(resultData))
#write.csv(resultData, file=paste0("RERconverge.continuous.fullEuth.GOenrichment.top250accOfAllAcc.translated.keepEarlyStops.noTraitWinsor.csv"), row.names=FALSE)
#write.csv(resultData, file=paste0("RERconverge.BINARY.noCOMP.lessStrict._GOenrichment.top100accOfAllAcc.translated.keepEarlyStops.noTraitWinsor.csv"), row.names=FALSE)
#write.csv(resultData, file=paste0("RERconverge.BINARY.noCOMP.lessStrict._GOenrichment.top50accOfAllAcc.translated.keepEarlyStops.noTraitWinsor.csv"), row.names=FALSE)

#Top accelerated out of all
total<-res$X[which(!(is.na(res$Rho)))]
inSubset<-c()
for(gene in total){
        if(gene %in% subset_acc){
                inSubset<-c(inSubset, 1)
        } else{
                inSubset<-c(inSubset, 0)
        }
}
if(length(subset_acc) != length(which(inSubset==1))){
        print("ERROR: length of 'subset_acc' file and 'inSubset' vector don't match up!")
}
names(inSubset)<-total
geneSelFunc<-function(iO){
        return(iO==1)
}
myGOdata<-new("topGOdata", description="Plt05Vtotal_acc", ontology="BP", allGenes=inSubset, geneSel=geneSelFunc, annot=annFUN.org, mapping="org.Hs.eg.db", ID="ensembl")
result<-runTest(myGOdata, algorithm = "classic", statistic = "fisher")
resultData<-GenTable(myGOdata, classFisher=result, orderBy="classFisher", ranksOf="classicFisher", topNodes=length(score(result)))
resultData$p.adj<-p.adjust(resultData$classFisher, method="BH")
print(head(resultData))
#write.csv(resultData, file=paste0("RERconverge.BINARY.noCOMP.lessStrict._GOenrichment.top250accOfAll.translated.keepEarlyStops.noTraitWinsor.csv"), row.names=FALSE)
#write.csv(resultData, file=paste0("RERconverge.BINARY.noCOMP.lessStrict._GOenrichment.top100accOfAll.translated.keepEarlyStops.noTraitWinsor.csv"), row.names=FALSE)
#write.csv(resultData, file=paste0("RERconverge.BINARY.noCOMP.lessStrict._GOenrichment.top50accOfAll.translated.keepEarlyStops.NoTraitWinsor.noAsk.csv"), row.names=FALSE)

#Decelerated only
total<-order_dec$X
inSubset<-c()
for(gene in total){
        if(gene %in% subset_dec){
                inSubset<-c(inSubset, 1)
        } else{
                inSubset<-c(inSubset, 0)
        }
}
if(length(subset_dec) != length(which(inSubset==1))){
        print("ERROR: length of 'subset_dec' file and 'inSubset' vector don't match up!")
}
names(inSubset)<-total
geneSelFunc<-function(iO){
        return(iO==1)
}
myGOdata<-new("topGOdata", description="Plt05Vtotal_dec", ontology="BP", allGenes=inSubset, geneSel=geneSelFunc, annot=annFUN.org, mapping="org.Hs.eg.db", ID="ensembl")
result<-runTest(myGOdata, algorithm = "classic", statistic = "fisher")
resultData<-GenTable(myGOdata, classFisher=result, orderBy="classFisher", ranksOf="classicFisher", topNodes=length(score(result)))
resultData$p.adj<-p.adjust(resultData$classFisher, method="BH")
print(head(resultData))
#write.csv(resultData, file=paste0("RERconverge.BINARY.noCOMP.lessStrict._GOenrichment.top250decOfAllDec.translated.keepEarlyStops.noTraitWinsor.csv"), row.names=FALSE)
#write.csv(resultData, file=paste0("RERconverge.BINARY.noCOMP.lessStrict._GOenrichment.top100decOfAllDec.translated.keepEarlyStops.noTraitWinsor.csv"), row.names=FALSE)
#write.csv(resultData, file=paste0("RERconverge.BINARY.noCOMP.lessStrict._GOenrichment.top50decOfAllDec.translated.keepEarlyStops.noTraitWinsor.noAsk.csv"), row.names=FALSE)


#Top decelerated out of all
total<-res$X[which(!(is.na(res$Rho)))]
inSubset<-c()
for(gene in total){
        if(gene %in% subset_dec){
                inSubset<-c(inSubset, 1)
        } else{
                inSubset<-c(inSubset, 0)
        }
}
if(length(subset_dec) != length(which(inSubset==1))){
        print("ERROR: length of 'subset_dec' file and 'inSubset' vector don't match up!")
}
names(inSubset)<-total
geneSelFunc<-function(iO){
        return(iO==1)
}
myGOdata<-new("topGOdata", description="Plt05Vtotal_dec", ontology="BP", allGenes=inSubset, geneSel=geneSelFunc, annot=annFUN.org, mapping="org.Hs.eg.db", ID="ensembl")
result<-runTest(myGOdata, algorithm = "classic", statistic = "fisher")
resultData<-GenTable(myGOdata, classFisher=result, orderBy="classFisher", ranksOf="classicFisher", topNodes=length(score(result)))
resultData$p.adj<-p.adjust(resultData$classFisher, method="BH")
print(head(resultData))
write.csv(resultData, file=paste0("RERconverge.continuous.rodentsubset.GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.csv"), row.names=FALSE)
#write.csv(resultData, file=paste0("RERconverge.BINARY.noCOMP.lessStrict._GOenrichment.top100decOfAll.translated.keepEarlyStops.noTraitWinsor.csv"), row.names=FALSE)
#write.csv(resultData, file=paste0("RERconverge.BINARY.noCOMP.lessStrict._GOenrichment.top50decOfAll.translated.keepEarlyStops.noTraitWinsor.noAsk.csv"), row.names=FALSE)

##Get GO spermatogenesis genes
#GO_spermatogenesis<-unlist(genesInTerm(myGOdata, "GO:0007283"))
#write(GO_spermatogenesis, file="GOspermatogenesis_genes.txt", ncol=1)
