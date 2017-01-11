#!/bin/bash
set -e

echo "${@:1}" > Input_Parameters.txt

#Inputs - Species 1
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
touch "PlantA_FA_"${FILENAMEA} 
touch "PlantA_GFF_"${FILENAMEB}



#Inputs - Species 2
# $7 - File, e.g. PlantB.fa
# $8 - File, e.g. PlantB.gff3
# $9 - String, e.g. ID=gene:
# $10 - Integer, e.g. 500/2000/5000
# $11 - --NoOverlap (if checked)
# $12 - --UseUTR (if checked)

fileC=$7
fileD=$8

FILENAMEC=${fileC##*/} # remove the path and leave only the file name
FILENAMED=${fileD##*/}

mv ${fileC} /apples/bin/utr_tool/
mv ${fileD} /apples/bin/utr_tool/

cd /apples/bin/utr_tool

./bo_utr.sh $7 $8 $9 $10 $11 $12

cp /apples/bin/utr_tool/works/promoters.fa /de-app-work/PlantB.fa
cp /apples/bin/utr_tool/works/promoters.bed /de-app-work/PlantB.bed
cp /apples/bin/utr_tool/works/utr3.bed /de-app-work/PlantB_utr3.bed
cp /apples/bin/utr_tool/works/utr5.bed /de-app-work/PlantB_utr5.bed


cd /de-app-work
# Make easier for the user to keep track of what input was used
touch "PlantB_FA_"${FILENAMEC} 
touch "PlantB_GFF_"${FILENAMED}




#Inputs - Conservation
# $13 - File, e.g. rbhSearch_result.txt
# $14 - Integer, window length, e.g. 30/[60]/80/100

if [ ! -d /apples/inputs ]; then mkdir /apples/inputs; fi
if [ ! -d /apples/outputs ]; then mkdir /apples/outputs; fi

fileE=$13

mv ${fileE} /apples/inputs/rbhSearch_result_PlantA_PlantB.txt
cp /de-app-work/* /apples/inputs/

cd /apples/bin

./conservationSearch_cyverse.pl -w $14
cp /apples/outputs/* /de-app-work/