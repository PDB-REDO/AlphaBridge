#!/bin/bash

set -e

trap 'error_handler $?' ERR

error_handler()
{
	echo "Error in find_interface: $1" | tee -a "${WORK}/find_interface.log"
	rm ${WORK}/running
	touch "${WORK}/error"
	exit 1
}

while getopts 'i:d:c:' OPTION; do
	case "$OPTION" in

		i)
			WORK="$OPTARG"
			;;
		d)
			DIR="$OPTARG"
			;;
		c)
			CONDA="$OPTARG"
			;;
		*)
			echo "usage: run-find-interface.sh -i <indir> -d <find_interface-dir> -c <mini-conda-dir>"
			exit 1
			;;
	esac
done

cd "${WORK}"

touch ${WORK}/running

(
	LANG=C eval "$(${CONDA}/bin/conda shell.bash hook)"
	LANG=C conda activate find_interface
	python3 "${DIR}/define_interfaces.py" -i "${WORK}"
) 2>&1 | tee -a "${WORK}/find_interface.log"

touch "${WORK}/finished"

rm ${WORK}/running
