#!/bin/bash

set -euo pipefail

help () {
	cat <<- EOF

	${0} [-kl] -i <input.bam> -f <forward.bed> -v <reverse.bed> -r <reference.fa>

	-i <file_id>	GDC file id.
	-s <case_id>	GDC case id.
	-k <token.tkn>	Path to GDC token file.
	
	-b <bam_name>	File name to be given to the downloaded GDC file.
	-d <output_dir>	Directory in which the results directory will be created. [Default=./]
	
	-f <forward.bed>	BED file of forward-coding regions in the capture kit.
	-v <reverse.bed>	BED file of reverse-coding regions in the capture kit.
	-r <reference.fa>	Reference genome FASTA file.

	-t <int>	Additional threads to assign. [Default=0]
	-q <base_quality>	Minimum base quality to be included in the pileup. [Default=37]

	-k	Do not count pileup bases that are a match. These counts are output as zeroes.

	EOF
}

threads=0
bq=37
out_dir="./"

while getopts ":i:s:d:b:r:f:v:t:k:q:kh" arg; do
	case "${arg}" in
		i)
			file_id="${OPTARG}"
			;;
		s)
			case_id="${OPTARG}"
			;;
		d)
			out_dir="${OPTARG}"
			;;
		b)
			bam_name="${OPTARG}"
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
		k)
			token="${OPTARG}"
			;;
		q)
			bq="${OPTARG}"
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

mkdir -p "${out_dir}/${case_id}/"
input_bam="${out_dir}/${case_id}/${bam_name}"
curl \
	--header "X-Auth-Token: $(cat "${token}")" \
	--output "${input_bam}" \
	"https://api.gdc.cancer.gov/slicing/view/${file_id}?region=chr21" 

if [[ "${input_bam}" == *bam ]]; then
	file_prefix="${input_bam%%.bam}"
elif [[ "${input_bam}" == *cram ]]; then
	file_prefix="${input_bam%%.cram}"
fi

filtered_bam="${file_prefix}.filtered.bam"

filter_wgs () {
	samtools view \
		--threads "${threads}" \
		--reference "${reference}" \
		--excl-flags 3852 \
		--min-MQ 60 \
		--remove-tag OQ \
		--input-fmt-option 'filter=ncigar==1 && ![XA] && ![SA]' \
		--bam \
		--output "${filtered_bam}" \
		"${input_bam}"

	samtools index "${filtered_bam}"
}

analyze_strand () {
	strand_name="${1}"
	bed="${2}"
	output_pileup="${filtered_bam%.bam}.${strand_name}.pileup"

	samtools mpileup \
		--fasta-ref "${reference}" \
		--no-BAQ \
		--min-BQ "${bq}" \
		--output-BP \
		--reverse-del \
		--positions "${bed}" \
		"${filtered_bam}" > "${output_pileup}"

	awk 'BEGIN {s=0} {s=s+$4} END {print s}' "${output_pileup}" > "${output_pileup%.pileup}.total_base_count.txt"

	awk '$5 !~ /^((\^.)?([\.,*#])(\$)?)+$/ {print}' "${output_pileup}" > "${output_pileup%.pileup}.no_match_positions.pileup"

	rm "${output_pileup}"
}

analyze_wgs () {
	output_pileup="${filtered_bam%.bam}.pileup"

	samtools mpileup \
		--fasta-ref "${reference}" \
		--no-BAQ \
		--min-BQ "${bq}" \
		--output-BP \
		--reverse-del \
		"${filtered_bam}" > "${output_pileup}"

	awk 'BEGIN {s=0} {s=s+$4} END {print s}' "${output_pileup}" > "${output_pileup%.pileup}.total_base_count.txt"

	awk '$5 !~ /^((\^.)?([\.,*#])(\$)?)+$/ {print}' "${output_pileup}" > "${output_pileup%.pileup}.no_match_positions.pileup"

	rm "${output_pileup}"

}

# Prefilter.
filter_wgs

# Do coding strand analysis.
(
	analyze_strand "forward" "${bed_forward}" &
	analyze_strand "reverse" "${bed_reverse}" &
	wait
	python pileup_parser.py \
		-id "${case_id}" \
		-pf "${file_prefix}.filtered.forward.no_match_positions.pileup" \
		-pr "${file_prefix}.filtered.reverse.no_match_positions.pileup" \
		-o "${file_prefix}.asym_table.tsv" \
		-s
) &

# Do WGS analysis.
(
	analyze_wgs
	python pileup_parser.py \
		-id "${case_id}" \
		-sp "${file_prefix}.filtered.no_match_positions.pileup" \
		-o "${file_prefix}.asym_table.wgs.tsv" \
		-s
) &

wait

rm "${input_bam}"* "${filtered_bam}"*
mkdir -p "${out_dir}/${case_id}/pileups/" "${out_dir}/${case_id}/total_base_counts/"
(cd "${out_dir}/${case_id}/" && tar -zcf "./pileups/${file_id}.pileups.tar.gz" ./*.no_match_positions.pileup)
rm "${out_dir}/${case_id}/"*.pileup
mv "${out_dir}/${case_id}/"*.total_base_count.txt "${out_dir}/${case_id}/total_base_counts/"

echo "Finished analyzing ${input_bam}"
