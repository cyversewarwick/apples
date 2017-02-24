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
		*)
			echo "Unkown option $1 with value $2"
			shift
		;;
	esac
	shift
done

echo "Species:    =${SPECIES}"
echo "Directory:  =${DBDIR}"
echo "Mode:       =${MODE}"
echo "Window:     =${WINDOW}"

ls -al ${DBDIR}

cd /apples/bin

if [ ! -d /apples/bin/tempfiles ]; then mkdir /apples/bin/tempfiles; fi

ulimit -c 0

# case ${MODE} in
# 	both)
# 		perl conservationSearch_cyverse_multiple.pl --species ${SPECIES} --dbdir ${DBDIR} -w ${WINDOW} || true
# 		perl conservationSearch_cyverse_multiple.pl --species ${SPECIES} --dbdir ${DBDIR} -w ${WINDOW} -p || true
# 	;;
# 	pseudo)
# 		perl conservationSearch_cyverse_multiple.pl --species ${SPECIES} --dbdir ${DBDIR} -w ${WINDOW} -p || true
# 	;;
# 	normal)
# 		perl conservationSearch_cyverse_multiple.pl --species ${SPECIES} --dbdir ${DBDIR} -w ${WINDOW} || true
# 	;;
# 	*)
# 		perl conservationSearch_cyverse_multiple.pl --species ${SPECIES} --dbdir ${DBDIR} -w ${WINDOW} || true
# 	;;
# esac

# cp /apples/outputs/* /de-app-work


echo "done"