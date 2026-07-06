#PURPOSE: Test if RTS is associated with mating system for each mammal clade

library(ape)
library(nlme)
library(phytools)
library(phylolm)
library(ggplot2)

clade<-"primate" #Options: rodent, primate, Laurasiatheria, all

print(paste("Testing for association between RTS and mating system for", clade))

#Read in Zoonomia tree
#Sourced from https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2.phast-242.nh
mytree<-read.tree("/ix1/nclark/shared_data/241mammals/241-mammalian-2020v2.phast-242.nh")

#Lists of mating systems that involve sperm competition or not
comp<-c("promiscuous", "polyandrous", "promiscuous; polyandrous", "polyandrous; promiscuous")
nocomp<-c("monogamous", "polygynous", "monogamous; polygynous", "polygynous; monogamous")

#Read in data and get clade
mydata<-read.csv("TableS1_MammalRTS_allPhenoData.csv", header=TRUE)
if(clade=="all"){
        cladedata<-mydata
} else{
        cladedata<-mydata[which(mydata$Major_clade==clade),]
}
print("Clade data info:")
print(dim(cladedata))
print(head(cladedata[,1:9]))

#Extract data
#If genetic mating system known, use that, otherwise use data based on social mating system/behavior
print("Running analysis using any mating system information...")
sperm_comp<-c() #A vector of 1s and 0s, where 1 means sperm competition and 0 means no sperm competition
rts<-c()
species<-c()
data_type<-c()
for(i in 1:nrow(cladedata)){
        if((!(is.na(cladedata$genetic_mating_system[i]))) && (cladedata$genetic_mating_system[i]!="")){
                if(cladedata$genetic_mating_system[i] %in% comp){
                        sperm_comp<-c(sperm_comp, 1)
                        species<-c(species, cladedata$Species_zoonomia[i])
                        rts<-c(rts, cladedata$rts[i])
                        data_type<-c(data_type, "genetic")
                } else if(cladedata$genetic_mating_system[i] %in% nocomp){
                        sperm_comp<-c(sperm_comp, 0)
                        species<-c(species, cladedata$Species_zoonomia[i])
                        rts<-c(rts, cladedata$rts[i])
                        data_type<-c(data_type, "genetic")
                } else{
                        print(paste("Genetic mating system is ambiguous for", cladedata$Species_zoonomia[i], "; skipping..."))
                }
        } else if(!(is.na(cladedata$mating_system[i]))){
                if(cladedata$mating_system[i] %in% comp){
                        sperm_comp<-c(sperm_comp, 1)
                        species<-c(species, cladedata$Species_zoonomia[i])
                        rts<-c(rts, cladedata$rts[i])
                        data_type<-c(data_type, "social")
                } else if(cladedata$mating_system[i] %in% nocomp){
                        sperm_comp<-c(sperm_comp, 0)
                        species<-c(species, cladedata$Species_zoonomia[i])
                        rts<-c(rts, cladedata$rts[i])
                        data_type<-c(data_type, "social")
		} else{
                        print(paste("Social mating system is ambiguous for", cladedata$Species_zoonomia[i], "; skipping..."))
                }
        } else{
                print(paste("No mating system data for", cladedata$Species_zoonomia[i], "; skipping..."))
        }
}

mydf<-as.data.frame(cbind(sperm_comp, rts, data_type))
rownames(mydf)<-species
mydf$rts<-as.numeric(as.character(rts))
print("Mating system data frame info:")
print(dim(mydf))
print(mydf)

#Run phylogenetic ANOVA
mydf_noNA<-mydf[which(!(is.na(mydf$rts))),]
pruned_tree<-keep.tip(mytree, rownames(mydf_noNA))
x<-mydf_noNA$sperm_comp
names(x)<-rownames(mydf_noNA)
y<-mydf_noNA$rts
names(y)<-rownames(mydf_noNA)
print(pruned_tree)
print(x)
print(y)
res<-phylANOVA(pruned_tree, x, y, p.adj="BH")
print(res)

#Plot as boxplot
pbox<-ggplot(mydf, aes(x=sperm_comp, y=rts)) + geom_boxplot(outlier.shape=NA) + theme_minimal()
pbox<-pbox + geom_jitter(shape=16, position=position_jitter(0.2))

pdf(paste("RTS_by_matingSystem.lessStrict",clade,"pdf", sep="."), onefile=TRUE)
print(pbox)
dev.off()

#Plot mating system data on tree
pdf(paste("matingSystemTree",clade,"pdf", sep="."), onefile=TRUE, height=20, width=10)
mating<-c()
for(i in mytree$tip.label){
        if(i %in% rownames(mydf)){
                mating<-c(mating, mydf$sperm_comp[which(rownames(mydf)==i)])
        } else{
                mating<-c(mating, NA)
        }
}
names(mating)<-mytree$tip.label
dotTree(mytree, mating, colors=setNames(c("blue","red","white"), c(0,1,NA)), fsize=0.5)
noNA<-mydf$sperm_comp
names(noNA)<-rownames(mydf)
dotTree(pruned_tree, noNA, colors=setNames(c("blue","red"), c(0,1)), fsize=0.5)
dev.off()
