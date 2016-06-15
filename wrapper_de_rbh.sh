#!/bin/bash
set -e

#inputs
# $1 - plantA.fa
# $2 - plantB.fa


if [ ! -d FASTADB ]; then mkdir FASTADB; fi
if [ ! -d BLASTDB ]; then mkdir BLASTDB; fi
if [ ! -d outputs ]; then mkdir outputs; fi

fileA=$1
fileB=$2

FILENAMEA=${fileA##*/} 
FILENAMEB=${fileB##*/}

mv ${fileA} /apples/FASTADB
mv ${fileB} /apples/FASTADB


cd /apples/bin
# The rbhSearch script accept both name.fa and name, but we remove the .fa here
perl rbhSearch.pl ${FILENAMEA%.fa} ${FILENAMEB%.fa} 

rm -r /apples/FASTADB/

rm -r /apples/BLASTDB/

cp /apples/outputs/* /de-app-work

echo "done"