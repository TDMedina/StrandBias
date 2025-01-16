#!/bin/bash

set -euo pipefail

help () {
	cat <<- EOF

	${0} -i <input.vcf.gz> -f <forward.bed> -v <reverse.bed>

	EOF
}


while getopts ":i:f:v:h" arg; do
	case "${arg}" in
		i)
			input_vcf="${OPTARG}"
			;;
		f)
			bed_forward="${OPTARG}"
			;;
		v)
			bed_reverse="${OPTARG}"
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

prepared="${input_vcf%.vcf.gz}.tumor_snvs.vcf.gz"
bcftools view --samples TUMOR --types snps --genotype ^miss --min-ac 1 -a "${input_vcf}" \
	| bcftools annotate --output-type z --output "${prepared}" --remove INFO/CSQ  
bcftools index "${prepared}"

forward_vcf="${prepared%.vcf.gz}.forward.vcf.gz"
reverse_vcf="${prepared%.vcf.gz}.reverse.vcf.gz"
bcftools view --regions-file "${bed_forward}" --output-type z --output "${forward_vcf}" "${prepared}"
bcftools view --regions-file "${bed_reverse}" --output-type z --output "${reverse_vcf}" "${prepared}"

echo "counts" > counts.tmp
{
	bcftools view -H "${forward_vcf}" | awk 'BEGIN {s=0} $4 == "G" && $5 == "T" {s=s+1} END {print s}'
	bcftools view -H "${forward_vcf}" | awk 'BEGIN {s=0} $4 == "C" && $5 == "A" {s=s+1} END {print s}'
	bcftools view -H "${reverse_vcf}" | awk 'BEGIN {s=0} $4 == "G" && $5 == "T" {s=s+1} END {print s}'
	bcftools view -H "${reverse_vcf}" | awk 'BEGIN {s=0} $4 == "C" && $5 == "A" {s=s+1} END {print s}'
} >> counts.tmp

cat << EOF > results.tmp
coding_region	mutation
forward	GtoT
forward	CtoA
reverse	GtoT
reverse	CtoA
EOF

paste results.tmp counts.tmp > counts.tsv
rm ./*.tmp ./*tumor_snvs*
