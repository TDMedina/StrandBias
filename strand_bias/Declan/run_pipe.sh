#!/bin/bash
# read one of the split files take first column

## Have python3 loaded
##sbatch run_pipe.sh ID1.txt ID1

while getopts ":i:n:" arg; do
	case "${arg}" in
		i)
			cram_list="${OPTARG}"
			;;
		n)
		 	id_dir="${OPTARG}"
		 	;;
		:)
		 	echo "Missing argument to option: ${OPTARG}"
		 	exit 1
		 	;;
		?)
		 	echo "Invalid option: ${OPTARG}"
		 	exit 1
		 	;;
	esac
done


be_polite () {
	while [[ $(squeue -u "${USER}" -states PD --noheader | wc -l) -gt 1000 ]]; do
		sleep 60
	done
}


while read -r case_id; do
	be_polite

	mkdir -p "pileups_bq37/${id_dir}/${case_id}"

	sbatch -J "${case_id}" -o "pileups_bq37/${id_dir}/${case_id}/${case_id}_slurmlogs_Tbases.txt" bases_mapped.sh "${case_id}" "${id_dir}"
done < "${cram_list}"
