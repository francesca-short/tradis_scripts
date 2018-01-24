# tradis_scripts
This repository contains some short utility scripts that I use in conjunction with the Bio-TraDIS pipeline.

tradis_insert_sites_FS.py
- Rewritten gene_insert_sites script from Bio-TraDIS pipeline. Main differences are that this script retains insert sites information for all feature classes except ones which you specifically exclude (lines 85-100), and there is a sanity check built in to ensure your plot file matches the length of the reference genome file you are using. Requires Biopython, gzip, numpy.

seq_saturation_test.py
- For examining the relationship between sequencing depth and number of insertion sites in your data. Takes a .bam and .bai file from your tradis mapping ("bacteria_tradis" command in the Bio-TraDIS pipeline) and calculates the number of unique insertion sites with different numbers of randomly samples reads. Requires pysam, matplotlib.
