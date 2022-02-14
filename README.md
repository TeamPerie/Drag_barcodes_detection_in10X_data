# lentiviral_barcode_detection_in10X_data


This is the code to extract and analyse lentiviral barcodes from 10x single cell RNA-seq data from the paper:
H3K27me3 is a determinant of chemotolerance in triple-negative breast cancer (Marsolier, Prompsy et al., Nature Genetics, 2022).

Originally the code was run on a cluster using the main script: run_10x_analysis_lenti_vallot.sh

If the code is to be run this way, there are several variables to be set at the beginning of the script (this is also all annotated in the script):

•	SCRIPT_DIR – directory where all the bash and R scripts are stored on the cluster
•	SAMPLES – sample names, separated by ;
•	NSAMP – number of samples
•	DATA_DIR – directory where separate 10x bam files are stored on cluster (each in a subdirectory SAMPLE_NAME/possorted_genome_bam.bam)
•	OUT_DIR – writable directory in which all output will go
•	R_EXEC – path to Rscript
•	Other parameters which can be changed, such as the number of bases to match in the flanking regions, which are currently set to the values used in the paper.

Alternatively, if you would like to the run the code locally, the variables can be set in the terminal for the bash parts of the code and the R scripts can be easily adapted to set the required variables at the start of each script.

For each sample, three scripts are run consecutively:
•	grep_flanking_10x_vbc2.sh – extracts reads which do not map to the transcriptome, matches a fixed sequence in the flanking regions of the lentiviral barcodes allowing one mismatch, and extracts the viral barcode, read ID, 10x cell barcode and 10x UMI. (Runs R codes agrep_flank_reads_left_vbc2.R and agrep_flank_reads_right_vbc2.R)
•	get_consensus_vbcs.sh – collates the viral barcodes which are extracted using matches to the up-/down-stream flanking regions in the same read, from different reads with the same UMI, and then from different UMIs from the same cell. A consensus approach is used, with the most common base returned for each position and then barcodes with poor consensus filtered. (Runs R code get_UMI_and_cell_consensus.R)
•	check_lenti_bcs_library.sh – cross-references barcodes output by previous step with the sequenced barcode library. Adds columns to output to say whether the detected viral barcode is in the library, and in how many reads of the two sequencing replicates of the library this barcode was found.


![image](https://user-images.githubusercontent.com/3630413/153838389-c783faff-a200-4675-a543-2c9df16ac82c.png)

