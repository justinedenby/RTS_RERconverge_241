#PURPOSE: create Newick trees with subsets of species -- in this example we made trees for each major clade using the breadTrees1() function 
library(RERconverge)
source("/ihome/nclark/emk270/MAMMAL_RTM_RER/breadtrees_original.R")


#read in our list of rodents 
rodents = read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/rodent_241_list.csv")
rodent_list = rodents$Species_zoonomia
print(rodent_list)

#make tree with our species subset
trees = breadTrees1("/ix1/nclark/ekopania/zoonomia_mammals_241/CODING_REGIONS/coding_region_RERtrees_filtered-50percent-20bp_translatedAA.keepEarlyStops_replaceAsterisk.tre", useSpecies=rodent_list)
#save tree as RDS
saveRDS(trees, file="/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/UPDATED_breadT_translated_rodents_translatedAA_keepEarlyStops_noAsk.rds")
print("saved rodent rds")


#read in our list of laurasiatherians
lauras = read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/laura_241_list.csv")
laura_list = lauras$Species_zoonomia
print(laura_list)

#make tree with our species subset
trees2 = breadTrees1("/ix1/nclark/ekopania/zoonomia_mammals_241/CODING_REGIONS/coding_region_RERtrees_filtered-50percent-20bp_translatedAA.keepEarlyStops_replaceAsterisk.tre", useSpecies=laura_list)

#save tree as RDS
saveRDS(trees2, file="/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/UPDATED_breadT_translated_lauras_translatedAA_keepEarlyStops_noAsk.rds")
print("saved laura rds")



#read in our list of primates
primates = read.csv("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/primate_241_list.csv")
primates_list = primates$Species_zoonomia
print(primates_list)

#make tree with our species subset
trees3 = breadTrees1("/ix1/nclark/ekopania/zoonomia_mammals_241/CODING_REGIONS/coding_region_RERtrees_filtered-50percent-20bp_translatedAA.keepEarlyStops_replaceAsterisk.tre", useSpecies=primates_list)

saveRDS(trees3, file="/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/UPDATED_breadT_translated_primates_translatedAA_keepEarlyStops_noAsk.rds")
print("saved primate rds")



