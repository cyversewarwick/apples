#!/bin/bash
set -e

echo "${@:1}" > Input_Parameters.txt

#Inputs

while [[ $# -gt 1 ]]; do
	key="$1"
	case $key in
		-s|--species)
			SPECIES="$2"
			shift
		;;
		-d|--dbdir)
			DBDIR="$2"
			shift
		;;
		-m|--mode)
			MODE="$2"
			shift
		;;
		-w|--window)
			WINDOW="$2"
			shift
		;;
		-a|--thresholda)
			THRESHOLDA="$2"
			shift
		;;
		-b|--thresholdb)
			THRESHOLDB="$2"
			shift
		;;
		*)
			echo "Unkown option $1 with value $2"
			shift
		;;
	esac
	shift
done

DEPATH=$(pwd)
DBFULLPATH=$(pwd)/${DBDIR}
OUTFULLPATH=$(pwd)/outputs

if [ ! -d ${OUTFULLPATH} ]; then mkdir ${OUTFULLPATH}; fi

echo "Species:    =${SPECIES}"
echo "Directory:  =${DBDIR}"
echo "Full Path:  =${DBFULLPATH}"
echo "Mode:       =${MODE}"
echo "Window:     =${WINDOW}"
echo "Thrshld_A:  =${THRESHOLDA}"
echo "Thrshld_B:  =${THRESHOLDB}"

ls -al ${DBFULLPATH}

ulimit -c 0

if [ ! -d /apples/bin/tempfiles ]; then mkdir /apples/bin/tempfiles; fi

cd /apples/bin

if [[ ${MODE} == 'normal' ]] || [[ ${MODE} == 'both' ]]; then
	echo "calling without p"
	perl conservationSearch_cyverse_multiple.pl \
		--species ${SPECIES} \
		--dbdir ${DBFULLPATH} \
		--outdir ${OUTFULLPATH} \
		--wsize ${WINDOW} \
		--tha ${THRESHOLDA} \
		--thb ${THRESHOLDB} || true
fi

if [[ ${MODE} == 'pseudo' ]] || [[ ${MODE} == 'both' ]]; then
	echo "calling with p"
	perl conservationSearch_cyverse_multiple.pl \
		--species ${SPECIES} \
		--dbdir ${DBFULLPATH} \
		--outdir ${OUTFULLPATH} \
		--wsize ${WINDOW} \
		--tha ${THRESHOLDA} \
		--thb ${THRESHOLDB} \
		--pseudo || true
fi

cd ${DEPATH}
tar -zcf con_${SPECIES//,/_}_w${WINDOW}_a${THRESHOLDA}_b${THRESHOLDB}_${MODE}.tar.gz outputs
# cp /apples/outputs/* /de-app-work

echo "done"