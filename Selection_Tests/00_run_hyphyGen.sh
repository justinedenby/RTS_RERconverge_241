#!/bin/bash
#SBATCH --job-name=hyphy_gen
#SBATCH --output=hyphy_gen-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jjd108@pitt.edu
#SBATCH --cluster=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200M

#Include hyphy scripts in path
export PATH=$PATH:/ix3/nclark/ekopania/MURINAE_REVISIONS/core/hyphy-interface/
echo $PATH

#PURPOSE: run hyphy_gen.py under our specific input sequences

#python ./hyphy-interface/hyphy_gen.py -i /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/spermato_0110filtered_noAsk-seq100-site100/nt -m relax -genetrees /ihome/nclark/jjd108/summer_2025/hyphy/early.stop.codon.fix_noAsk_PRUNED_GENE_TREES -o /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/relax/relax_earlystopfix_noAsk -tb /ihome/nclark/jjd108/summer_2025/for_hyphy_tenpercent_species.txt -p hyphy -part htc -tasks 64 -cpus 1 -mem 12000 -n relax_earlystopfix_noAsk

python ./hyphy-interface/hyphy_gen.py -i /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/spermato_0110filtered_noAsk-seq100-site100/nt -m busted-ph -genetrees /ihome/nclark/jjd108/summer_2025/hyphy/early.stop.codon.fix_noAsk_PRUNED_GENE_TREES -o /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/busted_ph/bustedph_earlystopfix_noAsk -tb /ihome/nclark/jjd108/summer_2025/for_hyphy_tenpercent_species.txt -p hyphy -part htc -tasks 64 -cpus 1 -mem 12000 -n bustedph_earlystopfix_noAsk



