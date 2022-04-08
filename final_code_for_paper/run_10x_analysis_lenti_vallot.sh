#!/bin/bash
#  run_10x_analysis_lenti_Vallot.sh
#
#  Created by A Lyne on 12/07/2020.
#  
# main script to run functions to extract lentiviral barcodes from 10x data
# 

#make id for analysis, used to name cluster jobs
ID="Vallot"
jid="_${ID}"

#TO DO: make directories as outlined below and insert paths
#any software that is called e.g. samtools in /bioinfo/local/build/Centos, needs to be added to PATH variable in .profile in your home directory on the cluster
# e.g. export PATH=$PATH:/bioinfo/local/build/Centos/samtools/samtools-1.9/bin

#ON CLUSTER
#MAKE DIRECTORY IN $HOME FOR SCRIPTS
#All scripts will be transferred here by scp
SCRIPT_DIR="/data/users/alyne/scrnaseq/lenti_bc_10x"

#names of samples, separated by ";"
SAMPLES="UNC_J33_6;5FU3_day202;5FU3_day50;5FU3_day77;5FU5_day171;5FU5_day67;5FU6_day214;5FU6_day33;5FU6_UNC_day33;DMSO3_day50;initial"
#number of samples
NSAMP=11

#MAKE BASE DIRECTORY FOR INPUT DATA (10x BAM FILES)
DATA_DIR="/data/tmp/scRNA_10x_lineage/"
mkdir $DATA_DIR

#IN DATA_DIR, MAKE A DIRECTORY FOR EACH SAMPLE (WITH THE NAME OF THE SAMPLE AS ABOVE)
#PLACE THE RELEVANT BAM FILE IN EACH DIRECTORY

#MAKE DIRECTORY FOR OUTPUT
OUT_DIR="/data/tmp/scRNA_10x_lineage_output"
mkdir "$OUT_DIR"


#ONCE ALL SCRIPTS HAVE BEEN TRANSFERRED TO SCRIPT_DIR ON CLUSTER
#CHANGE PERMISSIONS OF THIS FILE TO EXECUTABLE: chmod +x run_10x_analysis_lenti_vallot.sh

#cd to "$SCRIPT_DIR" and run using ./run_10x_analysis_lenti.sh

#once all qsubs have been run, transfer .csv output to local

#path to R version to call on cluster
R_EXEC="/bioinfo/local/build/Centos/R/R-3.4.0/bin/Rscript"


######################
## PARAMETERS TO SET
######################

#constant flanking sequences
FL_L="CGTCAACTAGAACACTCGAGATCAG"
FL_R="TGTGGTATGATGTATCATCTGAGTA"

#number of base pairs to match in flanking region of lentiviral barcode
#NB, max is 25 unless you add more bases to sequences in R agrep files
#shorter this is, the more barcodes but the more error...
NMATCH=20 #main match, one mismatch allowed
NMATCH_OPP=4 #exact match required on other side

#length of lentiviral barcode
BC_LENGTH=20

#flag, set to 1 to only output full length barcodes, 0 otherwise
FULL_LENGTH=1


#iterate through samples
for ((i=1;i<="$NSAMP";i++)); do

	#sample to analyse
	SAMPLE=$(echo "$SAMPLES" | awk -v l="$i" 'BEGIN {FS=";";OFS=";"} {print $l}')
	#location of bam file
	COUNT_DIR="${DATA_DIR}/${SAMPLE}"
	#directories for output
	AGREP_DIR="$OUT_DIR"/"$SAMPLE"/match_"$NMATCH"_"$NMATCH_OPP"
	mkdir "$OUT_DIR"/"$SAMPLE"
	mkdir "$AGREP_DIR"

	#get all unmapped reads to check for match to flanking region of barcode
	#grep for pattern match and output 10x barcode, UMI and 20bp viral barcode
	RUNID=$(qsub -l nodes=1:ppn=1 -l walltime=5:00:00 -l mem=5gb -N 10x_grep"$jid" -k oe -v SAMPLE="$SAMPLE",COUNT_DIR="$COUNT_DIR",NMATCH="$NMATCH",FL_L="$FL_L",FL_R="$FL_R",AGREP_DIR="$AGREP_DIR",R_EXEC="$R_EXEC",NMATCH_OPP="$NMATCH_OPP",BC_LENGTH="$BC_LENGTH",FULL_LENGTH="$FULL_LENGTH" grep_flanking_10x_vbc2.sh)


	#combine barcodes from left and right greps and output 10x cell barcode and UMI, plus viral barcodes
	RUNID2=$(qsub -W depend=afterany:"$RUNID" -l nodes=1:ppn=1 -l walltime=0:10:00 -l mem=1gb -N 10x_summ"$jid" -k oe -v SAMPLE="$SAMPLE",AGREP_DIR="$AGREP_DIR",R_EXEC="$R_EXEC",SCRIPT_DIR="$SCRIPT_DIR",NMATCH="$NMATCH",NMATCH_OPP="$NMATCH_OPP",BC_LENGTH="$BC_LENGTH" get_consensus_vbcs.sh)


done






















