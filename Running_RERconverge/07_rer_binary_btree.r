#PURPOSE: run RERconverge with binary trait AND subsetted Newick tree for each clade in our analysis (primate, rodents, laurasiatheria)
#requires helper function [preorder_fix.R] to preorder the phylogenetic tree to ensure the proper species alignment

library(RERconverge)
source("~/summer_2026/rerconverge_corrected_codon/preorder_fix.R")
 
#genetic mating system, non-COMP foreground, binary analysis
base_dir <- "/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon" 

#read in RDS file of clade specific Newick tree
#use patchTreesObjToPreorder() function to ensure that the tree is preordered
laura   <- readRDS(file.path(base_dir, "UPDATED_breadT_translated_lauras_translatedAA_keepEarlyStops_noAsk.rds"))
laura2  <- patchTreesObjToPreorder(laura)

#create object with species with no mating system related sperm competition to be used as the foreground
laura_nocomp_species <- c("Heterohyrax_brucei", "Procavia_capensis", "Bos_taurus",
                          "Hippopotamus_amphibius", "Hipposideros_galeritus",
                          "Miniopterus_schreibersii", "Scalopus_aquaticus",
                          "Sus_scrofa", "Ursus_maritimus", "Vulpes_lagopus")


#use RERconverge functions to create binary trait tree from foreground species
l_nocomp_tree <- foreground2Tree(laura_nocomp_species, laura2, clade = "all")
l_nocomp_phen <- tree2Paths(l_nocomp_tree, laura2)


#run RERconverge main function to create RER matrix with default parameters
l_mamRER <- RERconverge::getAllResiduals(laura2, transform = "sqrt")
#save RER matrix as RDS file
saveRDS(l_mamRER, file = file.path(base_dir, "l_mamRER_btree.rds"))

#run analysis as binary trait with default parameters
l_corrBin <- correlateWithBinaryPhenotype(l_mamRER, l_nocomp_phen,
                                          min.sp = 10, min.pos = 2, weighted = "auto")

#print the first 5 genes with the lowest P values
print(head(l_corrBin[order(l_corrBin$P), ]))

#write results to csv
write.csv(l_corrBin,
          file.path(base_dir, "LAURA_binary.genetic.mating.systems.nonCOMPforeground_results.csv"))





#the rest of the script follows the same detailing as above


# primate  <- readRDS(file.path(base_dir, "UPDATED_breadT_translated_primates_translatedAA_keepEarlyStops_noAsk.rds"))
# primate2 <- patchTreesObjToPreorder(primate)
# 
# primate_nocomp_species <- c("Aotus_nancymaae", "Callithrix_jacchus",
#                             "Gorilla_gorilla", "Nasalis_larvatus")
# 
# p_nocomp_tree <- foreground2Tree(primate_nocomp_species, primate2, clade = "all")
# p_nocomp_phen <- tree2Paths(p_nocomp_tree, primate2)
# 
# p_mamRER <- RERconverge::getAllResiduals(primate2, transform = "sqrt")
# saveRDS(p_mamRER, file = file.path(base_dir, "p_mamRER_btree.rds"))
# 
# p_corrBin <- correlateWithBinaryPhenotype(p_mamRER, p_nocomp_phen,
#                                           min.sp = 10, min.pos = 2, weighted = "auto")
# print(head(p_corrBin[order(p_corrBin$P), ]))
# write.csv(p_corrBin,
#           file.path(base_dir, "PRIMATE_binary.genetic.mating.systems.nonCOMPforeground_results.csv"))
# 
 
 

#rodent  <- readRDS(file.path(base_dir, "UPDATED_breadT_translated_rodents_translatedAA_keepEarlyStops_noAsk.rds"))
#rodent2 <- patchTreesObjToPreorder(rodent)

#rodent_nocomp_species <- c("Cavia_aperea", "Cavia_porcellus", "Cuniculus_paca",
#                           "Heterocephalus_glaber", "Perognathus_longimembris",
#                           "Mus_spretus", "Ondatra_zibethicus")

#r_nocomp_tree <- foreground2Tree(rodent_nocomp_species, x, clade = "all")
#r_nocomp_phen <- tree2Paths(r_nocomp_tree, x)

#r_mamRER <- RERconverge::getAllResiduals(rodent, transform = "sqrt")
#saveRDS(r_mamRER, file = file.path(base_dir, "r_mamRER_btree.rds"))

#r_corrBin <- correlateWithBinaryPhenotype(r_mamRER, r_nocomp_phen,
#                                          min.sp = 10, min.pos = 2, weighted = "auto")
#print(head(r_corrBin[order(r_corrBin$P), ]))
#write.csv(r_corrBin,
#          file.path(base_dir, "RODENT_binary.genetic.mating.systems.nonCOMPforeground_results.csv"))



 
 ##################################################
 
 
#euth  <- readRDS(file.path(base_dir, "coding_50percent-20bp_translated_keepEarlyStops_noAsk.rds"))
#euth_nocomp_species <- c("Heterohyrax_brucei", "Procavia_capensis", "Bos_taurus",
#                         "Hippopotamus_amphibius", "Hipposideros_galeritus",
#                         "Miniopterus_schreibersii", "Scalopus_aquaticus",
#                         "Sus_scrofa", "Ursus_maritimus", "Vulpes_lagopus","Aotus_nancymaae", "Callithrix_jacchus",
#                         "Gorilla_gorilla", "Nasalis_larvatus","Cavia_aperea", "Cavia_porcellus", "Cuniculus_paca",
#                         "Heterocephalus_glaber", "Perognathus_longimembris",
#                         "Mus_spretus", "Ondatra_zibethicus")

#e_nocomp_tree <- foreground2Tree(euth_nocomp_species, euth, clade = "all")
#e_nocomp_phen <- tree2Paths(e_nocomp_tree, euth)

#euth_species = read.csv(file.path(base_dir, "TableS1_MammalRTS_allPhenoData - Sheet1.csv"))
#e_mamRER <- RERconverge::getAllResiduals(euth, transform = "sqrt", useSpecies =euth_species$Species_zoonomia)
#saveRDS(e_mamRER, file = file.path(base_dir, "updated_euth_mamRER.rds"))

#e_corrBin <- correlateWithBinaryPhenotype(e_mamRER, e_nocomp_phen,
#                                          min.sp = 10, min.pos = 2, weighted = "auto")
#print(head(e_corrBin[order(e_corrBin$P), ]))
#write.csv(e_corrBin,
#          file.path(base_dir, "euth_binary.genetic.mating.systems.nonCOMPforeground_results.csv"))
