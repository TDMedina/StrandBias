#!/bin/bash

set -euo pipefail

help () {
	cat <<- EOF

	${0} -m <manifest> -f <forward.bed> -v <reverse.bed>

	EOF
}


while getopts ":m:f:v:h" arg; do
	case "${arg}" in
		m)
			manifest="${OPTARG}"
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

counter=1
counter_total="$(cat "${input_vcf_list}" | wc -l )"

cat << EOF > counts.tsv
forward	forward	forward	forward	reverse	reverse	reverse	reverse
GT	GT	CA	CA	GT	GT	CA	CA
TOTAL	PASS	TOTAL	PASS	TOTAL	PASS	TOTAL	PASS
EOF

echo ""

while read -r -a input_vcf; do

	project_id="${input_vcf[6]}"
	case_id="${input_vcf[5]}"
	file_id="${input_vcf[1]}"
	file_name="${input_vcf[2]}"
	file_path="./samples/${case_id}/${file_id}/${file_name}"

	echo -ne "\r${counter}/${counter_total}"
	prepared="${file_path}%.vcf.gz}.tumor_snvs.vcf.gz"
	
	bcftools view \
	--samples TUMOR \
	--types snps \
	--genotype ^miss \
	--min-ac 1 \
	--trim-alt-alleles \
	--min-alleles 2 \
	--max-alleles 2 \
	"${file_path}" \
	| bcftools annotate --output-type z --output "${prepared}" --remove INFO/CSQ
	bcftools index "${prepared}"

	forward_vcf="${prepared%.vcf.gz}.forward.vcf.gz"
	reverse_vcf="${prepared%.vcf.gz}.reverse.vcf.gz"
	bcftools view --regions-file "${bed_forward}" --output-type z --output "${forward_vcf}" "${prepared}"
	bcftools view --regions-file "${bed_reverse}" --output-type z --output "${reverse_vcf}" "${prepared}"

	declare -A counts
	counts["forward_gt"]=$(bcftools view -H "${forward_vcf}" | awk 'BEGIN {s=0} $4 == "G" && $5 == "T" {s=s+1} END {print s}')
	counts["forward_gt_pass"]=$(bcftools view -H -f PASS "${forward_vcf}" | awk 'BEGIN {s=0} $4 == "G" && $5 == "T" {s=s+1} END {print s}')
	counts["forward_ca"]=$(bcftools view -H "${forward_vcf}" | awk 'BEGIN {s=0} $4 == "C" && $5 == "A" {s=s+1} END {print s}')
	counts["forward_ca_pass"]=$(bcftools view -H -f PASS "${forward_vcf}" | awk 'BEGIN {s=0} $4 == "C" && $5 == "A" {s=s+1} END {print s}')
	counts["reverse_gt"]=$(bcftools view -H "${reverse_vcf}" | awk 'BEGIN {s=0} $4 == "G" && $5 == "T" {s=s+1} END {print s}')
	counts["reverse_gt_pass"]=$(bcftools view -H -f PASS "${reverse_vcf}" | awk 'BEGIN {s=0} $4 == "G" && $5 == "T" {s=s+1} END {print s}')
	counts["reverse_ca"]=$(bcftools view -H "${reverse_vcf}" | awk 'BEGIN {s=0} $4 == "C" && $5 == "A" {s=s+1} END {print s}')
	counts["reverse_ca_pass"]=$(bcftools view -H -f PASS "${reverse_vcf}" | awk 'BEGIN {s=0} $4 == "C" && $5 == "A" {s=s+1} END {print s}')

	echo -e "${counts["forward_gt"]}\t${counts["forward_gt_pass"]}\t${counts["forward_ca"]}\t${counts["forward_ca_pass"]}\t${counts["reverse_gt"]}\t${counts["reverse_gt_pass"]}\t${counts["reverse_ca"]}\t${counts["reverse_ca_pass"]}" >> counts.tsv
	((counter++))

done < "${input_vcf_list}"
