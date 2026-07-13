#PURPOSE: run RERconverge with a continous trait with trees that have subsets of species -- we specifically subsetted the major clades present in our analysis.  
library(RERconverge)

#RDS files are the Newick trees created from the breadTree source function
laura = readRDS("UPDATED_breadT_translated_lauras_translatedAA_keepEarlyStops_noAsk.rds")
primate = readRDS("UPDATED_breadT_translated_primates_translatedAA_keepEarlyStops_noAsk.rds")
rodent = readRDS("UPDATED_breadT_translated_rodents_translatedAA_keepEarlyStops_noAsk.rds")

#read in data for analysis
trusted_data = read.csv("/ihome/nclark/jjd108/semester2/teste_trait_finding/log_rts_csv_trusted_241continuous.csv")

#run RERconverge main function with default parameters to get RER matrix for each clade, save as RDS to eliminate recalculating
	#bread_RER = getAllResiduals(laura, transform = "sqrt", weighted = T, scale = T)
	#bread_RER = getAllResiduals(rodent, transform = "sqrt", weighted = T, scale = T)
bread_RER = getAllResiduals(primate, transform = "sqrt", weighted = T, scale = T)
saveRDS(bread_RER, file = "/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/RER_breadT_TRANSLATED_primates_KeepEarlyStops_noAsk.rds")

#create dataframe for input data
z241_df = data.frame(trusted_data$Species_zoonomia, trusted_data$LOG_rts)
z241_vector <- setNames(z241_df$trusted_data.LOG_rts, z241_df$trusted_data.Species_zoonomia)

# continuous data to input tree 
charP = char2Paths(z241_vector, primate)
print(charP)

#runs continuous trait analysis with clade specific RER matrix
cphen = correlateWithContinuousPhenotype(bread_RER, charP, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)


#Rho, p-value, and Benjamini-Hochberg corrected p-value
cphen$stat = -log10(cphen$P) * sign(cphen$Rho)
head(cphen[order(cphen$stat, decreasing = TRUE),])

#save RERconverge results as a csv file -- repeated for all clades (rodents, primates, laurasiatherians) 
genOutput = write.csv(cphen, file = "./translated_breadTree_primates_only_KeepEarlyStops_noAsk_RERresults.csv")

