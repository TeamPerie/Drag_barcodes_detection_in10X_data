#!/bin/bash
#  split_all_reads_10x.sh
#
#  Created by A Lyne on 26/04/2018.
#

cd "$AGREP_DIR"

#get unmapped reads from 10x bam file
samtools view -f 4 "$COUNT_DIR"/possorted_genome_bam.bam > all_reads_"${SAMPLE}".txt


#get sequences to grep from flanking regions
NCHAR=${#FL_L}
START=$((NCHAR-NMATCH+1))
GREP_L=$(echo $FL_L | cut -c "$START"-"$NCHAR")

NCHAR=${#FL_R}
GREP_R=$(echo $FL_R | cut -c 1-"$NMATCH")


#make regex to match exactly or with one mismatch
patt1=$GREP_L
for ((i=0; i<${#patt1[0]}; i++)); do 
patt1+=( "${patt1[0]:0:i}.${patt1[0]:i+1}" )
done
regex1=$(IFS='|'; echo "${patt1[*]}")

patt2=$GREP_R
for ((i=0; i<${#patt2[0]}; i++)); do 
patt2+=( "${patt2[0]:0:i}.${patt2[0]:i+1}" )
done
regex2=$(IFS='|'; echo "${patt2[*]}")

regex=${regex1}"|"${regex2}


#take reads which either match left or right flanking region
awk -v regex="$regex" '$0 ~ regex' all_reads_"${SAMPLE}".txt > reads_pattern_match_"${SAMPLE}".txt


# take reads which have10x BC and UMI, sequence, read name,
awk 'BEGIN {FS="\t";OFS="\t"} { for (i=10; i<=NF; ++i) { for (j=10; j<=NF; ++j) { if ($i ~ "CB:" && $j ~ "UB:") print substr($i,6),$10,$1,substr($j,6) } } }' reads_pattern_match_"${SAMPLE}".txt > unaligned_10xbc_seq_"$SAMPLE".txt


#run R code to grep NMATCH bp on one side and NMATCH_OPP on the other side of barcode
R_FILES="$PBS_O_WORKDIR/agrep_flank_reads_right_vbc2.R"
WD="$AGREP_DIR"
INFILE=unaligned_10xbc_seq_"$SAMPLE".txt

#grep on right, 3' end of barcode
OUTFILE="$AGREP_DIR"/agrep_10xbc_and_vbc_"${SAMPLE}"_"$NMATCH"_"$NMATCH_OPP"_right.txt
$R_EXEC $R_FILES $WD $INFILE $OUTFILE $BC_LENGTH $NMATCH $NMATCH_OPP $GREP_R $GREP_L $FULL_LENGTH

#grep on left, 5' end of barcode
R_FILES="$PBS_O_WORKDIR/agrep_flank_reads_left_vbc2.R"
OUTFILE="$AGREP_DIR"/agrep_10xbc_and_vbc_"${SAMPLE}"_"$NMATCH"_"$NMATCH_OPP"_left.txt
$R_EXEC $R_FILES $WD $INFILE $OUTFILE $BC_LENGTH $NMATCH $NMATCH_OPP $GREP_R $GREP_L $FULL_LENGTH


#uncomment if you want to keep this output
# rm all_reads_"${SAMPLE}".txt reads_pattern_match_"${SAMPLE}".txt unaligned_10xbc_seq_"$SAMPLE".txt unaligned_10xbc_seq_"$SAMPLE"_cells.txt










