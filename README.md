# tradis_scripts
This repository contains some short utility scripts that I use in conjunction with the Bio-TraDIS pipeline.

tradis_insert_sites_FS.py
- Rewritten gene_insert_sites script from Bio-TraDIS pipeline. Main differences are that this script retains insert sites information for all feature classes except ones which you specifically exclude (lines 85-100), and there is a sanity check built in to ensure your plot file matches the length of the reference genome file you are using. Requires Biopython, gzip, numpy.

seq_saturation_test.py
- For examining the relationship between sequencing depth and number of insertion sites in your data. Takes a .bam and .bai file from your tradis mapping ("bacteria_tradis" command in the Bio-TraDIS pipeline) and calculates the number of unique insertion sites with different numbers of randomly sampled reads. Requires pysam, matplotlib.

find_2direction_tn5.py
- Tn5 transposons can flip in place once inserted, which shows up in the plot file as two insertion sites on opposite strands, 9bp apart (this is due to target site duplication during the transposition). This script finds the proportion of these sites in your plot file, and collapses reads to a single site.

tradis_insert_steps_py3.py
- similar to tradis_insert_sites_FS.py, but assigns reads and insertion sites based on genome location (in user-defined chunks) rather than features. Useful for making charts or finding interesting intergenic regions.

Tradis_comparison_positive_selection.R
- Modified version of tradis_comparison.R script from BioTraDIS pipeline, with additional filtering by insertion index per gene in hit calling. This is for experiments with strong positive selection such as phage or >MIC antibiotic treatment, where calling hits based on read count comparisons alone may give spurious hits due to secondary mutations unrelated to the transposon insertion (when a gene shows increased read counts at one or two insertion sites, but loss of other insertion sites, this is presumed to be due to secondary mutations). Runs in R, requires getopt, edgeR, dplyr.
