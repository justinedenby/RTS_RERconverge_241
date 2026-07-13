#PURPOSE: run RERconverge with RTS as a continuous trait

library(RERconverge)

#read in input data for continuous analysis
teste_data = read.csv("/ihome/nclark/jjd108/semester2/teste_trait_finding/log_rts_csv_trusted_241continuous.csv")
print(teste_data)
print("read in csv")

#input species for analysis -- subset of 62/241
species = readLines("/ihome/nclark/jjd108/semester2/teste_trait_finding/241_specieslist.txt")

print("verify species: ")
print(species)


#first instance of running analysis -- must run readTrees, then save the tree as an RDS file to load faster 
#tree = readTrees("/ix1/nclark/ekopania/zoonomia_mammals_241/CODING_REGIONS/coding_region_RERtrees_filtered-50percent-20bp_translatedAA.keepEarlyStops_replaceAsterisk.tre")
#print("read translated tree")
#saveRDS(tree, file = "/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/coding_50percent-20bp_translated_keepEarlyStops_noAsk.rds")
tree = readRDS("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/coding_50percent-20bp_translated_keepEarlyStops_noAsk.rds")



#main RERconverge function, also save as RDS 
mamRER = getAllResiduals(tree, useSpecies=species, transform = "sqrt", weighted = T, scale = T)
saveRDS(mamRER, file = "/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/rer_res_matrix_corrected_earlystop_noAsk.rds")

#make input data as a proper vector object
teste_content = teste_data$LOG_rts
names(teste_content) = teste_data$Species_zoonomia
print(teste_content)


#char2Paths function correlates the RER matrix with your continous data -- RTS in this context
charP = char2Paths(teste_content, tree)


#correlation with continuous phenotype with default parameters
res = correlateWithContinuousPhenotype(mamRER, charP, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)


#RERconverge results output for continuous test: csv contains Gene, Rho, p-value, and Benjamini-Hochberg corrected p-values for analysis

genOutput = write.csv(res, file = "/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/RTS_RER_translated_keepEarlyStops_noAsk_results.csv")

print("order by P")
head(res[order(res$P),])

