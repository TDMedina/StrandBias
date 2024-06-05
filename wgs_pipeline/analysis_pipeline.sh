#!/bin/bash

set -euo pipefail

help () {
	cat <<- EOF

	${0} [-kl] -i <input.bam> -f <forward.bed> -v <reverse.bed> -r <reference.fa>

	-i <input.bam>	BAM file to analyze.
	-f <forward.bed>	BED file of forward-coding regions in the capture kit.
	-v <reverse.bed>	BED file of reverse-coding regions in the capture kit.
	-r <reference.fa>	Reference genome FASTA file.

	-t <int>	Additional threads to assign. [Default=0]

	-k	Do not count pileup bases that are a match. These counts are output as zeroes.
	-l	Parse pileups and write temporary tables sequentially to save memory. No summary
		tables are output in this mode.

	EOF
}

low_memory=""
skip_match_bases=""
threads=0
bq=37

while getopts ":s:i:r:f:v:t:q:klh" arg; do
	case "${arg}" in
		s)
			sample_id="${OPTARG}"
			;;
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
		t)
			threads="${OPTARG}"
			;;
		q)
			bq="${OPTARG}"
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


if [[ "${input_bam}" == *bam ]]; then
	file_prefix="${input_bam%%.bam}"
elif [[ "${input_bam}" == *cram ]]; then
	file_prefix="${input_bam%%.cram}"
fi

analyze_strand () {
	strand_name="${1}"
	bed="${2}"
	output_bam="${file_prefix}.filtered.${strand_name}.bam"

	samtools view \
	--threads "${threads}" \
	--reference "${reference}" \
	--region-file "${bed}" \
	--threads "${threads}" \
	--reference "${reference}" \
	--excl-flags 3852 \
	--min-MQ 60 \
	--remove-tag OQ \
	--input-fmt-option 'filter=ncigar==1 && ![XA] && ![SA]' \
	--bam \
	--output "${output_bam}" \
	"${input_bam}"

	output_pileup="${file_prefix}.filtered.${strand_name}.pileup"
	samtools mpileup \
	--fasta-ref "${reference}" \
	--no-BAQ \
	--min-BQ "${bq}" \
	--output-BP \
	--reverse-del \
	"${output_bam}" > "${output_pileup}" 2> "${output_pileup}.stderr"

	awk 'BEGIN {s=0} {s=s+$4} END {print s}' "${output_pileup}" > "${output_pileup%.pileup}.total_base_count.txt"

	awk '$5 !~ /^((\^.)?([\.,*#])(\$)?)+$/ {print}' "${output_pileup}" > "${output_pileup%.pileup}.no_match_positions.pileup"

	rm "${output_pileup}"
}

analyze_strand "forward" "${bed_forward}" &
analyze_strand "reverse" "${bed_reverse}" &

wait

python pileup_parser.py \
	-id "${sample_id}" \
	-pf "${file_prefix}.filtered.forward.no_match_positions.pileup" \
	-pr "${file_prefix}.filtered.reverse.no_match_positions.pileup" \
	-o "${file_prefix}.asym_table.tsv" \
	-s \
	${skip_match_bases}


# Clean up.
dest_dir="$(dirname "${input_bam}")"
file_id="$(basename "${input_bam}")"
file_id="${file_id%%.bam}"

bash 05.cleanup.sh -d "${dest_dir}" -f "${file_id}"
rm "${input_bam}"
