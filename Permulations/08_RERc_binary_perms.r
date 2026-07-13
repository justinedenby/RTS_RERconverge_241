#PURPOSE: run permulations on RERconverge binary trait analysis
library(RERconverge)

rts_data<-read.csv("~/summer_2026/rerconverge_corrected_codon/enrichment/mammal_rtsAndMatingSystem_zoonomia241mammals.withGeneticMating_high.ish_confidence_and.logrts.csv.", header=TRUE)


#define the foreground species
nocomp_species = c("Procavia_capensis","Bos_taurus","Camelus_dromedarius","Canis_lupus","Felis_catus","Hippopotamus_amphibius","Hipposideros_galeritus","Miniopterus_schreibersii",
                   "Odobenus_rosmarus","Vulpes_lagopus","Scalopus_aquaticus","Sus_scrofa",
                   "Aotus_nancymaae","Callithrix_jacchus","Gorilla_gorilla","Castor_canadensis","Cavia_aperea","Cavia_porcellus","Cuniculus_paca","Fukomys_damarensis",
                   "Heterohyrax_brucei","Perognathus_longimembris","Microtus_ochrogaster",
                   "Mus_spretus","Ondatra_zibethicus")

#read in Newick tree and RERconverge matrix from getAllResiduals()
tree = readRDS("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/coding_50percent-20bp_translated_keepEarlyStops_noAsk.rds")
mamRER = readRDS("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/rer_res_matrix_corrected_earlystop_noAsk.rds")

#create proper paths between foreground and binary trait
nocomp_tree = foreground2Tree(nocomp_species, tree, clade="all")
nocomp_phen  = tree2Paths(nocomp_tree, tree)   # <-- this is the numeric paths vector

#run binary analysis
corrBin = correlateWithBinaryPhenotype(mamRER, nocomp_phen, min.sp = 10, min.pos = 2,
        				weighted = "auto")

head(corrBin[order(corrBin$P),])
corrBin$stat = -log10(corrBin$P) * sign(corrBin$Rho)
head(corrBin[order(corrBin$stat, decreasing = TRUE),])



#binary permulations necessities -- denote which species are sister species on tree; only really necessary for enrichment analyses
sisters = list("clade1"= c("Procavia_capensis","Heterohyrax_brucei"),
	       "clade2"=c("Microtus_ochrogaster","Ondatra_zibethicus"))
root_sp = "Homo_sapien"
masterTree = tree$masterTree

#run 5000 permulations
permCC = getPermsBinary(5000, nocomp_species, sisters, root_sp, mamRER, tree,
			masterTree, permmode="cc",calculateenrich = F)

permpvalCC = permpvalcor(corrBin, permCC)

head(permpvalCC)
colnames(permpvalCC)

#align to corrBin by gene name and pull the numeric p-value
corrBin$permpval    = permpvalCC[match(rownames(corrBin), rownames(permpvalCC)), "permpval"]
corrBin$permpvaladj = p.adjust(corrBin$permpval, method = "BH")

head(corrBin[order(corrBin$permpval), ])

#save permulation results to csv
write.csv(corrBin,
          file = "/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/WITH_5000_PERMS_binary_fulleuth_RTS_RER_translated_keepEarlyStops_noAsk_results.csv")
