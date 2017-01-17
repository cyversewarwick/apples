#!/bin/bash
set -e

echo "${@:1}" > Input_Parameters.txt

#inputs
# $1 - plantA.fa
# $2 - plantB.fa
# $3 - simple/forked


if [ ! -d /apples/FASTADB ]; then mkdir /apples/FASTADB; fi
if [ ! -d /apples/BLASTDB ]; then mkdir /apples/BLASTDB; fi
if [ ! -d /apples/outputs ]; then mkdir /apples/outputs; fi

fileA=$1
fileB=$2

FILENAMEA=${fileA##*/} # remove the path and leave only the file name
FILENAMEB=${fileB##*/}

mv ${fileA} /apples/FASTADB/PlantA.fa
mv ${fileB} /apples/FASTADB/PlantB.fa

cd /apples/bin

# The rbhSearch script accept both name.fa and name, but we remove the .fa here
# ${FILENAMEB%.fa} <- this was supposed to remove the ".fa" suffix but sometimes this could be .fasta or anything else
# So we ditch this filter (next line of comment) and use a more generic species name to make sure everything runs.
#perl rbhSearch.pl ${FILENAMEA%.fa} ${FILENAMEB%.fa}
# perl rbhSearch.pl PlantA PlantB

case ${3} in
	simple)
	perl rbhSearch.pl PlantA PlantB
	;;
	forked)
	if [ ! -d /apples/tempdir ]; then mkdir /apples/tempdir; fi
	perl rbhSearchForked.pl PlantA PlantB
	;;
	*)
	perl rbhSearch.pl PlantA PlantB
	;;
esac

# rm -r /apples/FASTADB/

# rm -r /apples/BLASTDB/

cp /apples/outputs/* /de-app-work

cd /de-app-work
touch "PlantA_"${FILENAMEA} # So that the user can keep track of what input was used
touch "PlantB_"${FILENAMEB}

echo $(nproc)
echo "done"