#PURPOSE: run RERconverge permulations on continuous trait RERconverge runs
library(RERconverge)

#Read in CSV
teste_csv = read.csv("/ihome/nclark/jjd108/semester2/teste_trait_finding/log_rts_csv_trusted_241continuous.csv")
teste_content = teste_csv$LOG_rts
names(teste_content) = teste_csv$Species_zoonomia
print(teste_content)


#read in Newick Tree and RERconverge matrix from getAllResiduals()
tree = readRDS("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/coding_50percent-20bp_translated_keepEarlyStops_noAsk.rds")
mamRER = readRDS("/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/rer_res_matrix_corrected_earlystop_noAsk.rds")

charP = char2Paths(teste_content, tree)
res = correlateWithContinuousPhenotype(mamRER, charP, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)

#midpoint root tree for permulations
mt = tree$masterTree
r_mt = midpoint_root(mt)


#run permulations without enrichment
permsnoenrich = RERconverge::getPermsContinuous(5000, teste_content, mamRER, NULL,
                                                tree, r_mt, calculateenrich = F)
corpermpvals  = RERconverge::permpvalcor(res, permsnoenrich)
 
##attach permulation p-values to RERconverge matrix, aligned by GENE NAME
if (is.null(dim(corpermpvals))) {
  res$permpval = corpermpvals[match(rownames(res), names(corpermpvals))]
} else {
  res$permpval = corpermpvals[match(rownames(res), rownames(corpermpvals)), 1]
}
 
## keep raw and adjusted p-values in separate columns
res$permpvaladj = p.adjust(res$permpval, method = "BH")
 
## save results to csv
write.csv(res, file = "/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/WITH_5000_PERMS_RTS_RER_translated_keepEarlyStops_noAsk_results.csv")
