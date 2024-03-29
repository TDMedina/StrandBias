#!/bin/bash

set -euo pipefail

help () {
	cat <<- EOF

	${0} -i <input.bam> -f <forward.bed> -v <reverse.bed> -r <reference>

	EOF
}

while getopts ":i:r:f:v:h" arg; do
	case "${arg}" in
		i)
			input_bam="${OPTARG}"
			;;
		f)
			bed_forward="${OPTARG}"
			;;
		v)
			bed_reverse="${OPTARG}"
			;;
		r)
			reference="${OPTARG}"
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

split_by_read_orientation () {
	local orientations=("F1" "R2" "F2" "R1")
	for ori in "${orientations[@]}"; do
		local ori_subset="${2%%.bam}.${ori}.bam"
		{ bash 03.split_by_read_orientation.sh -d "${ori}" -o "${ori_subset}" -r "${1}" -i "${2}" \
			&& { bash 04.pileup.sh -r "${reference}" -i "${ori_subset}" & \
				bash 00.flagstat.sh -i "${ori_subset}" & } } &
	done
}

filtered="${input_bam%%.bam}.filtered.bam"

# Baseline flagstat in background.
bash 00.flagstat.sh -i "${input_bam}" &

# 1: Filter BAM, then flagstat result in background.
bash 01.filter_bam.sh -t 5 -i "${input_bam}" -r "${reference}" -o "${filtered}" \
	&& { 00.flagstat.sh -i "${filtered}" & }

# 2-4: Subset by strand, subset by orientation, and run pileup.
declare -A strands
strands["forward"]="${bed_forward}"
strands["reverse"]="${bed_reverse}"

for strand in "${!strands[@]}"; do
	region_subset="${filtered%%.bam}.${strand}_coding.bam"
	{ bash 02.subset_bam_by_region.sh \
		-i "${filtered}" \
		-r "${reference}" \
		-b "${strands[${strand}]}" \
		-o "${region_subset}" \
		&& { bash 00.flagstat.sh -i "${region_subset}" & \
			split_by_read_orientation "${reference}" "${region_subset}" & } } &
done
wait

sleep 10

for strand in "${!strands[@]}"; do
	orientations=("F1R2" "F2R1")
	for ori in "${orientations[@]}"; do
		merged="${filtered%.bam}.${strand}_coding.${ori}.bam"
		bash 03b.merge_read_orientations.sh \
			-r "${reference}" \
			-a "${filtered%.bam}.${strand}_coding.${ori:2}.bam" \
			-b "${filtered%.bam}.${strand}_coding.${ori::2}.bam" \
			-o "${merged}" \
			&& bash 04.pileup.sh -r "${reference}" -i "${merged}" -o "${merged%.bam}.pileup" &
	done
done
wait

sleep 10

file_prefix=$(readlink -f "${input_bam}")
file_prefix="${file_prefix%%.bam}"
python pileup_parser.py \
	-f "${file_prefix}" \
	-o "${input_bam%%.bam}.asym_table.tsv" \
	--summary-path "${input_bam%%.bam}.asym_summary_table.tsv" \
	-s \
	-m "combined"

pileup_dir="$(dirname $(readlink -f ${input_bam}))/pileups/"
mkdir -p "${pileup_dir}"
mv "${input_bam%%.bam}"*.no_match_positions.pileup "${pileup_dir}"
rm "${file_prefix}"*bam "${input_bam%%.bam}"*pileup "${input_bam%%.bam}"*bai
