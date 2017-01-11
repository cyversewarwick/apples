#!/bin/bash
set -e

# Example ./bo_utr.sh inputs/Arabidopsis_thaliana.TAIR10.31.dna_rm.toplevel.fa inputs/Arabidopsis_thaliana.TAIR10.31.gff3 ID=gene: 500 --No --UseUTR

#inputs
# $1 - File, e.g. PlantA.fa
# $2 - File, e.g. PlantA.gff3
# $3 - String, e.g. ID=gene:
# $4 - Integer, e.g. 500/2000/5000
# $5 - --NoOverlap (if checked)
# $6 - --UseUTR (if checked)

fileA=$1
fileB=$2

FILENAMEA=${fileA##*/} # remove the path and leave only the file name
FILENAMEB=${fileB##*/}

mv ${fileA} /apples/bin/utr_tool/
mv ${fileB} /apples/bin/utr_tool/

cd /apples/bin/utr_tool

./bo_utr.sh $1 $2 $3 $4 $5 $6

cp /apples/bin/utr_tool/works/promoters.fa /de-app-work/PlantA.fa
cp /apples/bin/utr_tool/works/promoters.bed /de-app-work/PlantA.bed
cp /apples/bin/utr_tool/works/utr3.bed /de-app-work/PlantA_utr3.bed
cp /apples/bin/utr_tool/works/utr5.bed /de-app-work/PlantA_utr5.bed


cd /de-app-work
# Make easier for the user to keep track of what input was used
touch "FA_"${FILENAMEA} 
touch "GFF_"${FILENAMEB}
echo "${@:1}" > Input_Parameters.txt
