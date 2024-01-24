#!/bin/bash

set -euo pipefail

help () {
	cat <<- EOF

	${0} [-kl] -i <input.bam> -f <forward.bed> -v <reverse.bed> -r <reference.fa>

	-i <input.bam>	BAM file to analyze.
	-f <forward.bed>	BED file of forward-coding regions in the capture kit.
	-v <reverse.bed>	BED file of reverse-coding regions in the capture kit.
	-r <reference.fa>	Reference genome FASTA file.

	-k	Do not count pileup bases that are a match. These counts are output as zeroes.
	-l	Parse pileups and write temporary tables sequentially to save memory. No summary
		tables are output in this mode.

	EOF
}

low_memory=""
skip_match_bases=""

while getopts ":i:r:f:v:klh" arg; do
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
		l)
			low_memory="--low-memory"
			;;
		k) 
			skip_match_bases="--skip-match-bases"
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
	&& { bash 00.flagstat.sh -i "${filtered}" & }

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

# Merge F and R reads to make F1R2 and F2R1 bams.
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

# Parse pileups to make summary tables of F1R2 and F2R1 counts.
# file_prefix=$(readlink -f "${input_bam}")
# file_prefix="${file_prefix%%.bam}"
file_prefix="${input_bam%%.bam}"


python pileup_parser.py \
	-f "${file_prefix}" \
	-o "${file_prefix}.asym_table.tsv" \
	--summary-path "${file_prefix}.asym_summary_table.tsv" \
	-s \
	-m "combined" \
	${skip_match_bases} \
	${low_memory}


# Clean up.
dest_dir="$(dirname "${input_bam}")"
file_id="$(basename "${input_bam}")"
file_id="${file_id%%.bam}"

bash 05.cleanup.sh -d "${dest_dir}" -f "${file_id}"

# mkdir -p "${dest_dir}/flagstats/" "${dest_dir}/logs/" "${dest_dir}/pileups/" "${dest_dir}/total_base_counts/"
# mv "${dest_dir}/"*.flagstat "${dest_dir}/flagstats/"
# mv "${dest_dir}/"*.std* "${dest_dir}/logs/"
# tar -zcf "${dest_dir}/pileups/${file_id}.pileups.tar.gz" "${dest_dir}/"*.pileup
# rm "${dest_dir}/"*.pileup
# mv "${dest_dir}/"*.total_base_count.txt "${dest_dir}/total_base_counts/"

# rm "${file_prefix}"*.bam "${input_bam%%.bam}"*.bai
