#!/bin/bash
#PURPOSE: Use newick utils (nw_prune) to prune gene trees to match species in alignment AFTER removing species with early stop codons
#
# Job name:
#SBATCH --job-name=prune_gene_trees
#SBATCH --output=prune_gene_tree-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jjd108@pitt.edu
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=2G #Not sure if I should mess with these...
#
## Command(s) to run:
module load python/ondemand-jupyter-python3.11
source activate /ihome/nclark/emk270/software/envs/hyphy2.5.62

#CHANGE THIS - output directory for your pruned gene trees
mkdir "early.stop.codon.fix_noAsk_PRUNED_GENE_TREES"
#CHANGE THIS - directory with post-filtering alingment files (i.e., fastas with early stop codon sequences removed)
ls /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/spermato_0110filtered_noAsk-seq100-site100/nt | while read file; do
	#Get just the gene name from the full file path and filename
	gene=$(echo ${file} | cut -d "/" -f -1 | cut -d "." -f 1)
	echo "Pruning tree for ${gene}..."
	#Get the names of species present in the post-filtering alignment fasta and writes them to a file
	#CHANGE THIS - file path should match the "ls" statement above (path to iqtree output fastas)
	grep ">" /ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/spermato_0110filtered_noAsk-seq100-site100/nt/${file} > temp_keep_tips.txt
	sed -i 's/> //' temp_keep_tips.txt
	#CHANGE THIS PATH - get the IQtree gene tree for this gene
	tree="/ihome/nclark/jjd108/summer_2026/hyphy_post_translation_corrections/hyphy_iqtree_gen_noAsk/loci/${gene}/${gene}.treefile"
	#Make a new subdirectory for this gene where the pruned treefile will be output
	mkdir early.stop.codon.fix_noAsk_PRUNED_GENE_TREES/${gene}
	#Prunes all tips EXCEPT the ones present in the given file
	#-v makes it so we are KEEPING the tips in temp_keep_tips.txt
	#-f is the treefile
	nw_prune -v -f ${tree} temp_keep_tips.txt > "early.stop.codon.fix_noAsk_PRUNED_GENE_TREES/${gene}/${gene}.treefile"
done

echo "Done!"
