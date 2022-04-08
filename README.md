# Lentiviral barcode extraction from 10x scRNA-seq data
--------------

This is the code to extract and analyse lentiviral barcodes from 10x single cell RNA-seq data as used in the paper:

H3K27me3 is a determinant of chemotolerance in triple-negative breast cancer (Marsolier, Prompsy et al., Nature Genetics, 2022).

Full details of the approach are given in the methods section of the paper.


### Running the code

Originally the code was run on a cluster using the main script: run_10x_analysis_lenti_vallot.sh

If the code is to be run this way, there are several variables to be set at the beginning of the script (this is also annotated in the script):

* SCRIPT_DIR – directory where all the bash and R scripts are stored on the cluster
* SAMPLES – sample names, separated by ;
* NSAMP – number of samples
* DATA_DIR – directory where separate 10x bam files are stored on cluster (each in a subdirectory SAMPLE_NAME/possorted_genome_bam.bam)
* OUT_DIR – writable directory in which all output will go
* R_EXEC – path to Rscript
* Other parameters which can be changed, such as the number of bases to match in the flanking regions, which are currently set to the values used in the paper.

Alternatively, if you would like to the run the code locally, the variables can be set in the terminal for the bash parts of the code and the R scripts can be easily adapted to set the required variables at the start of each script.


### Requirements

* Samtools
* R (libraries seqinr, Biostrings)


### How does it work?

For each sample, two scripts are run consecutively:
* grep_flanking_10x_vbc2.sh – extracts reads which do not map to the transcriptome, matches a fixed sequence in the flanking regions of the lentiviral barcodes allowing one mismatch, and extracts the viral barcode, read ID, 10x cell barcode and 10x UMI. (Runs R codes agrep_flank_reads_left_vbc2.R and agrep_flank_reads_right_vbc2.R)
* get_consensus_vbcs.sh – collates the viral barcodes which are extracted using matches to the up-/down-stream flanking regions in the same read, from different reads with the same UMI, and then from different UMIs from the same cell. A consensus approach is used, with the most common base returned for each position and then barcodes with poor consensus filtered. (Runs R code get_UMI_and_cell_consensus.R)


### Input and example output

ScRNA-seq data from Marsolier, Prompsy et al (2022) can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164715). The main inputs to the lentiviral barcode extraction pipeline are the 10x bam files with the aligned and non-aligned reads. These can be found by clicking on the GSM link for a given sample, then the SRX link and finally the SRR link. You can find example output for sample "initial" in file consensus_cell_10xbc_and_vbc_initial_20_4.csv.zip above.


