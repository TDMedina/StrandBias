#!/bin/bash


set -euo pipefail

help () {
	cat <<- EOF

	${0} -d <directory> -f <file_id>

	EOF
}

while getopts ":d:r:f:v:h" arg; do
	case "${arg}" in
		d)
			dest_dir="${OPTARG}"
			;;
		f)
			file_id="${OPTARG}"
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

mkdir -p "${dest_dir}/flagstats/" "${dest_dir}/logs/" "${dest_dir}/pileups/" "${dest_dir}/total_base_counts/"
mv "${dest_dir}/"*.flagstat "${dest_dir}/flagstats/"
mv "${dest_dir}/"*.std* "${dest_dir}/logs/"
(cd "${dest_dir}" && tar -zcf "./pileups/${file_id}.pileups.tar.gz" ./*.pileup)
rm "${dest_dir}/"*.pileup
mv "${dest_dir}/"*.total_base_count.txt "${dest_dir}/total_base_counts/"
rm "${dest_dir}/${file_id}"*.bam "${dest_dir}/${file_id}"*.bai
