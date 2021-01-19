#!/bin/bash
#  produce_output_summary.sh
#
#  Created by A Lyne on 18/06/2020.
#
# script to get consensus viral barcode for each UMI and to combine greps from left and right

cd "$AGREP_DIR"

#put all output into one file and tag with whether full pattern match was on right or left of barcode
awk 'BEGIN {FS=OFS="\t"} {print $0,"L"}' agrep_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP"_left.txt > agrep_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP".txt
awk 'BEGIN {FS=OFS="\t"} {print $0,"R"}' agrep_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP"_right.txt >> agrep_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP".txt


#order by cell barcode then UMI
sort -k1,2 agrep_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP".txt > agrep_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP"_sorted.txt


#R script to get consensus viral barcode for each UMI and then for each cell
R_FILE="$SCRIPT_DIR"/get_UMI_and_cell_consensus.R

$R_EXEC $R_FILE $AGREP_DIR agrep_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP"_sorted.txt consensus_UMI_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP".csv consensus_cell_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP".csv $BC_LENGTH


#remove uneccessary files
# rm agrep_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP"_left.txt agrep_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP"_right.txt agrep_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP".txt