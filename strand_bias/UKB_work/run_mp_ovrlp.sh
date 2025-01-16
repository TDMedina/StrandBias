#!/bin/bash

set -uo pipefail

###~~~~~
# There are multiple parts to this command firstly I am calling pileup. I am forced to used the deprecated samtools version as pileup has been migrated to bcftools and no longers has the functionality I require. I am imposing as basq of 20, mapq of 60 & max depth of 200. -x keeps the overlapping reads -O outputs the position in reads and --ouput-Qname out puts read ids.
#
# The output is then piped to a depth filter {20} and then piped through a indel filter. Before overlapping bases are called via script python
#

# awk filters depth and loci that contain indels and loci that contain non reference.
# the sed filter on the single read pileups_bq37 removes read positional information and the read mapping quality succeeding the ^ (^! is 33- 33= 0,^] is 93-33). If you know your decimal representation of the ascii table this makes sense. This filter essentially forces the match string to be equal to the quality string.

1>&2 echo "Running sample ${1}"
mkdir -p "pileups_bq37/${2}/${1}"

# get single and overlapping read pileups_bq37
#~/bin/samtools/samtools mpileup --input-fmt-option 'filter=![XA] && ![SA]' -F 3852 -B --ff SUPPLEMENTARY,DUP,UNMAP,SECONDARY,QCFAIL -Q 37 -q 60 -d 200 -f /tmp/GRCh38_full_analysis_set_plus_decoy_hla.fa -x --output-QNAME -O /tmp/${1}.cram | tee >(awk '$4 > 20 && $5!~/\+|\-[0-9]+[ACGTNacgtn]+/ && $5 ~ /[a-zA-Z]+/' | python parse_mup3.py pileups_bq37/${2}/${1}/${1}.test.gz pileups_bq37/${2}/${1}) >(awk '$4>20 && $5!~/\+|\-[0-9]+[ACGTNacgtn]+/ && $5~/[ACGTNacgtn]+/' | sed 's/\^//g; s/\$//g; s/\]//g; s/\!//g; s/<//g' | awk '$4 == length($5) {print $1,$2,$3,$4,$5,$7}' | bgzip -c > pileups_bq37/${2}/${1}/${1}.test.mpup.gz) >/dev/null

~/bin/samtools/samtools mpileup \
	--input-fmt-option 'filter=ncigar==1 && ![XA] && ![SA]' \
	-F 3852 \
	-B \
	--ff SUPPLEMENTARY,DUP,UNMAP,SECONDARY,QCFAIL \
	-Q 37 \
	-q 60 \
	-d 200 \
	-f /tmp/GRCh38_full_analysis_set_plus_decoy_hla.fa \
	-x \
	--output-QNAME \
	-O "/tmp/${1}.cram" \
	| awk '$4 > 20 && $5!~/\+|\-[0-9]+[ACGTNacgtn]+/ && $5 ~ /[a-zA-Z]+/' \
	| python parse_mup3.py "pileups_bq37/${2}/${1}/${1}.test.gz" "pileups_bq37/${2}/${1}"

if ! [ $? -eq 0 ]; then
	echo "mpileup overlap failed" >> "pileups_bq37/${2}/${1}/${1}.error"
	exit 1
fi

~/bin/samtools/samtools mpileup \
	--input-fmt-option 'filter=ncigar==1 && ![XA] && ![SA]' \
	-F 3852 \
	-B \
	--ff SUPPLEMENTARY,DUP,UNMAP,SECONDARY,QCFAIL \
	-Q 37 \
	-q 60 \
	-d 200 \
	-f /tmp/GRCh38_full_analysis_set_plus_decoy_hla.fa \
	-O "/tmp/${1}.cram" \
	| awk '$4>20 && $5!~/\+|\-[0-9]+[ACGTNacgtn]+/ && $5~/[ACGTNacgtn]+/' \
	| sed 's/\^//g; s/\$//g; s/\]//g; s/\!//g; s/<//g' \
	| awk '$4 == length($5) {print}' \
	| bgzip -c > "pileups_bq37/${2}/${1}/${1}.mpup.gz"

if ! [ $? -eq 0 ]; then
	echo "mpileup failed" >> "pileups_bq37/${2}/${1}/${1}.error"
	exit 1
fi


#FYI 3852 can be decoded here https://broadinstitute.github.io/picard/explain-flags.html

## Remaining analyses
## get total mapped reads. might be better to count the exome capture regions
~/bin/samtools/samtools view \
	-c \
	-q 60 \
	-f 3 \
	-F 3852 \
	--input-fmt-option 'filter=![XA] && ![SA]' \
	-T /tmp/GRCh38_full_analysis_set_plus_decoy_hla.fa \
	"/tmp/${1}.cram" > "pileups_bq37/${2}/${1}/${1}_MAPPED_cnt.txt"

if ! [ $? -eq 0 ]; then
	echo "read count failed" >> "pileups_bq37/${2}/${1}/${1}.error"
	exit 1
fi


# count supplementary reads, i.e XA tag as the bitwise operator doesnt have it. Supplementary reads are reads that also elsewhere
count=$(~/bin/samtools/samtools view \
	-c \
	-F 3852 \
	--input-fmt-option 'filter=[XA]' \
	-T /tmp/GRCh38_full_analysis_set_plus_decoy_hla.fa \
	"/tmp/${1}.cram")

if ! [ $? -eq 0 ]; then
	echo "Supplemental failed" >> "pileups_bq37/${1}/${1}.error"
	exit 1
fi

echo "${1} ${count}" > "pileups_bq37/${2}/${1}/${1}_Supplementary_reads_cnt.txt"


# chimeric reads write to out- These are secondary or read fragments where the remaining linear sequence on the read maps else where. process CRAM later
~/bin/samtools/samtools view \
	-F 3852 \
	--input-fmt-option 'filter=[SA]' \
	-T /tmp/GRCh38_full_analysis_set_plus_decoy_hla.fa \
	-C \
	-o "pileups_bq37/${2}/${1}/${1}_chimeras.cram" \
	"/tmp/${1}.cram"

if ! [ $? -eq 0 ]; then
	echo "Chimeras failed" >> "pileups_bq37/${2}/${1}/${1}.error"
	exit 1
fi


# count reads mapping to exome capture region.
bam2bed < <(~/bin/samtools/samtools view \
		-F 3852 \
		--input-fmt-option 'filter=![XA] && ![SA]' \
		-h \
		-q 60 \
		-T /tmp/GRCh38_full_analysis_set_plus_decoy_hla.fa /tmp/${1}.cram) \
	| bedmap \
		--echo \
		--count \
		--skip-unmapped \
		--fraction-map 0.01 \
		/tmp/Exome-IDT_V1.bed -\
	| bgzip -c > "pileups_bq37/${2}/${1}/${1}_Exome_capture_regions_cnt.txt.gz"

if ! [ $? -eq 0 ]; then
	echo "Exome capture failed" >> "pileups_bq37/${2}/${1}/${1}.error"
	exit 1
fi

# insert Brians overlapping fragment code here

~/bin/samtools/msCaller2 "/tmp/${1}.cram" "/data3/declan/scratch/ms_shortlist.bed" "/tmp/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
 | sed 's/\\//g' \
 | bgzip -c >  "pileups_bq37/${2}/${1}/${1}_mscaller.txt.gz"

if ! [ $? -eq 0 ]; then
	echo "MS overlap failed" >> pileups_bq37/${1}/${1}.error
	exit 1
fi
