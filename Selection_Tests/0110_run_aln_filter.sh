#!/bin/bash
#SBATCH --job-name=aln_filter
#SBATCH --output=aln_filter-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jjd108@pitt.edu
#SBATCH --cluster=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-02:00:00 #Time limit 2hours
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB

module load biopython/1.73
module load python/bioconda-3.7-2019.03
module load python/ondemand-jupyter-python3.11
source activate /ihome/nclark/emk270/software/envs/python3.11
echo "source activated"

#Include Gregg's core scripts in path
#CHANGE THIS TO YOUR CWD WHICH SHOULD ALSO CONTAIN THE corelib DIRECTORY
export PATH=$PATH:/ihome/nclark/jjd108/summer_2025/hyphy/corelib/
echo $PATH

#Run aln_filter.py
#-i is the directory of input fasta alignments (nucleotide sequences)
#-s 100 sets gappy filter to 100% (i.e., up to 100% gaps tolerated, no sequences will be removed due to gaps)
#-c 100 sets codon window threshold to 100% (i.e., no sequences will be removed due to codon filter)
#-o is the output directory (the script will create this directory, so it should not exist yet when you run the script)
#python aln_filter.py -i /ihome/nclark/jjd108/summer_2025/hyphy/confirmed_spermato_codon/ -s 100 -c 100 -o /ihome/nclark/jjd108/summer_2025/hyphy/spermato_confirmed_codon/
#python aln_filter.py -i /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/confirmed_spermato_corrected_codon/ -s 100 -c 100 -o /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/spermato_filtered/
python aln_filter.py -i /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/spermato_corrected_codon_No_Ask -s 100 -c 100 -o /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/spermato_0110filtered_noAsk

