#!/bin/bash


set -euo pipefail

help () {
	cat <<- EOF

	${0} -i <file_id> -n <file_name> -s <case_id>

	-i <file_id>	GDC file id.
	-n <file_name>	GDC file name.
	-s <case_id>	GDC case id.

	EOF
}

while getopts ":i:n:s:h" arg; do
	case "${arg}" in
		i)
			file_id="${OPTARG}"
			;;
		n)
			file_name="${OPTARG}"
			;;
		s)
			case_id="${OPTARG}"
			;;
		h)
			help
			exit 0
			;;
		:)
			echo "Missing required argument for option: -${OPTARG}"
			help
			exit 1
			;;
		?)
			echo "Invalid option: -${OPTARG}"
			help
			exit 1
			;;
	esac
done

if [ ${OPTIND} -eq 1 ]; then help; exit 0; fi

asym_table="./TestSample/${case_id}/${file_name%.bam}.asym_table.wgs.tsv"

if [ -f "${asym_table}" ]; then
	echo "WGS table already exists."
	exit 0
fi

pileup_dir="./TestSample/${case_id}/pileups/"
pileup_tar="${pileup_dir}/${file_id}.pileups.tar.gz"
pileup_file="./${file_name%.bam}.filtered.no_match_positions.pileup"

tar --directory "${pileup_dir}" -x "${pileup_file}" -f "${pileup_tar}"
python pileup_parser.py \
	-id "${case_id}" \
	-sp "${pileup_dir}/${pileup_file}" \
	-o "${asym_table}" \
	-s \
	--skip-match-bases

rm "${pileup_dir}/${pileup_file}"
echo "WGS table created for case ID: ${case_id}"
