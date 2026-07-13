#!/bin/bash
#PURPOSE: Run Rscript for toGO analyses
#
# Job name:
#SBATCH --job-name=topGO
#SBATCH --output=topGO-%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cluster=htc
##SBATCH --partition=<partition>
##SBATCH --constraint=amd,genoa
#SBATCH --mail-user=jjd108@pitt.edu
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00 #Time limit 1 days
#SBATCH --qos=normal
#SBATCH --mem-per-cpu=4G
#
## Command(s) to run:

module load gcc/12.2.0 r/4.4.0

#Rscript 00b_append_geneIDs.r

#Runs topGO
Rscript 01_mammal_testes_GO.r

#Rscript 02_primate_testes_GO.r

#Rscript 04_compare_murine_primate_GO.r

#Generate enrichment density plots
#Rscript 05_enrichmentDensityPlots.r
#Rscript 05_enrichmentDensityPlots.byStat.r

#Rscript 06_mammal_testes_GO.permPs.r

#Rscript 07_compare_taxa.PortsAndMurines.r

#Looks for overlaps and correlations among taxonomic groups
#Rscript 07_compare_taxa.r

#module load python/ondemand-jupyter-python3.11
#source activate /ihome/nclark/jjd108/SETUP/envs/my_env
#python 03_clean_go_csvs.py

#Rscript 08_matingSystem_RTS_corr.r

echo "Done!"

