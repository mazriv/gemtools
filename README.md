


This project contains scripts to run tools from the gemtools suite to support Brain TRI project.

The original gemtools github and further info can be found at:
https://github.com/sgreer77/gemtools . (citation is also below)

Cancer genome variants: useful for single nucleotide variants (SNVs), copy number variants (CNVs), structural variants (SVs), fusions, haplotyping. 

Topics: linked-reads whole genome sequencing (LR-WGS), whole genome sequencing (WGS), barcodes, gems, 10x Genomics.

Type: scripts.





Gemtools contains the following tools to analyze 10X Genomics linked-reads whole genome sequencing data:

get_phased_basic: Obtain phasing information for all SNVs in the vcf file

get_phase_blocks: Summarize phase blocks -- coordinates, size, number of phased heterozygous SNVs per phase block etc.

get_phased_bcs: For a particular phase block, return the haplotype 1 and haplotype 2 barcodes

get_bcs_in_region: Get all the barcodes that exist in a given region(s) of the genome

count_bcs_list: Determine presence and quantity of given barcodes across a given region

plot_hmw: Generate a plot of the mapping locations of reads with each barcode

plot_vars_and_blocks: For a particular region, plot the heterozygous variants and phase blocks

set_bc_window: Generate windows around SV breakpoints for SV analysis

get_shared_bcs: Determine barcodes shared between SV breakpoints

set_hap_window: Generate windows around SV breakpoints for haplotype analysis

assign_sv_haps: Assign SV barcodes to existing haplotypes (SNVs)

count_bcs: Determine presence and quantity of given barcodes across a given region surrounding the SV breakpoints

plot_hmw: Generate a plot of the mapping locations of reads with each barcode (SAME AS ABOVE)

extract_reads: Obtain reads with particular barcodes from Long Ranger fastq files (where fastq output is R1,R2,I1)

extract_reads_interleaved: Obtain reads with particular barcodes from Long Ranger fastq files (where fastq output is RA,I1,I2)



Citation for gemtools can be found at: https://academic.oup.com/bioinformatics/article/35/21/4397/5426055





