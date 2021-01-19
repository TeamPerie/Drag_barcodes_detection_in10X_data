#!/bin/bash
#  split_all_reads_10x.sh
#
#  Created by A Lyne on 26/04/2018.
#

cd "$AGREP_DIR"

#check if barcodes are in filtered library
awk 'NR == FNR { id[$1]; next } {if($2 in id){OFS="\t"; print $0,"1"}else{OFS="\t"; print $0,"0"}}' <(awk 'NR%2==0' ../../additional_files/LG22_filtered.fa) <(awk 'BEGIN {FS=",";OFS="\t"} {$1=$1}1' consensus_cell_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP".csv) | awk 'NR>1' > in_lib.txt


#check how many reads detected in replicate 1 and replicate 2 of sequencing
awk 'NR == FNR { a[$1] = $1; b[$1] = $3; next } {if($2 in a){OFS="\t"; print $0,b[$2]}else{OFS="\t"; print $0,"0"}}' <(awk 'BEGIN {FS=OFS="\t"} NR>1' ../../additional_files/A890R1_R1-res-all.txt) in_lib.txt > in_seq_rep1.txt

awk 'NR == FNR { a[$1] = $1; b[$1] = $3; next } {if($2 in a){OFS="\t"; print $0,b[$2]}else{OFS="\t"; print $0,"0"}}' <(awk 'BEGIN {FS=OFS="\t"} NR>1' ../../additional_files/A890R2_R1-res-all.txt) in_seq_rep1.txt > in_seq_rep2.txt


#add flag to say whether observed barcodes are in lentiviral library and how many reads they are observed with in two sequencing replicates of barcode library
echo "BC_10x,cons_BC_lenti,Numi_per_base,Prop_match_cons,side,kept,in_lib,seq_rep1,seq_rep2" > consensus_cell_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP"_annot.csv
awk 'BEGIN {FS="\t";OFS=","} {$1=$1}1' in_seq_rep2.txt >> consensus_cell_10xbc_and_vbc_"$SAMPLE"_"$NMATCH"_"$NMATCH_OPP"_annot.csv


rm in_lib.txt in_seq_rep1.txt in_seq_rep2.txt


