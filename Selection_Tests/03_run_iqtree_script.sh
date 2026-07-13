#!/bin/bash
##PURPOSE: Run Gregg's iqtree_gen.py script to generate files and scripts for making iqtree concat tree
#
# Job name:
#SBATCH --job-name=iqtree_gen
#SBATCH --output=iqtree_gen.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jjd108@pitt.edu # Where to send mail
#SBATCH --cluster=htc
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem=1M
#
## Command(s) to run:

#-i Path to input alignments (amino acids, filtered)
#-o Path to output directory
#-n Name for this run
#Other parameters should be fine as is, but python 04_iqtree_gen.py -h will tell you what they all are, as well as other options
#python 04_iqtree_gen.py -tasks 64 -i /ihome/nclark/jjd108/summer_2025/hyphy/spermato_genes -o ./iqtree_TEST -p iqtree2 -n iqtree_TEST -part htc
#python 04_iqtree_gen.py -i /ihome/nclark/jjd108/semester4/ift88_fel/ -o /ihome/nclark/jjd108/summer_2026/TEST_iqtree_ift88 -p iqtree2 -n TEST_iqtree_ift88 -part htc
#python 04_iqtree_gen.py -tasks 64 -i /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/confirmed_spermato_corrected_codon -o /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/hyphy_iqtree_gen -p iqtree -n hyphy_iqtree_gen -part htc

python 04_iqtree_gen.py -tasks 64 -i /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/spermato_corrected_codon_No_Ask -o /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/hyphy_iqtree_gen_noAsk -p iqtree -n hyphy_iqtree_gen_noAsk -part htc



echo "Done!" 
