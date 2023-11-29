#!/bin/bash
#SBATCH --partition=normal,highmem,MSC
#SBATCH --job-name=filter_bam


#SBATCH --time=00-00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

#SBATCH --error=%x.%j.stderr
#SBATCH --output=%x.%j.stdout

#SBATCH --mail-user=tyler.medina@universityofgalway.ie
#SBATCH --mail-type=BEGIN,END,FAIL

# Start Script

function help {
	cat <<- EOF

	01.filter_bam.sh [-r <reference.fa> -o <output.bam> -p] -i <input.bam>

	-i <input.bam>	Input BAM file.
	-o <output.bam>	Output BAM file. [Default=<input>.filtered.bam]
	-r <reference.fasta>	Reference FASTA file. [Default=./reference.fa]
	-t <int>	Additional threads to assign. [Default=0]
	-p	If set, script will echo the output name for pipeline purposes.

	EOF
}

reference="./reference.fa"
pipe_echo=0
threads=0

while getopts ":r:i:o:t:ph" arg; do
	case "${arg}" in
		h)
			help
			exit 0
			;;
		i)
			input_bam="${OPTARG}"
			;;
		o)
			output_bam="${OPTARG}"
			;;
		r)
			reference="${OPTARG}"
			;;
		t)
			threads="${OPTARG}"
			;;
		p)
			pipe_echo=1
			;;
		:)
			echo "Missing argument to option: -${OPTARG}"
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
if [ -z ${input_bam+x} ]; then echo "Missing required argument: -i"; help; exit 1; fi 

if [ -z ${output_bam+x} ]; then output_bam="${input_bam%%.bam}.filtered.bam"; fi 

samtools view \
	--threads "${threads}" \
	--reference "${reference}" \
	--bam \
	--output "${output_bam}" \
	--excl-flags 3852 \
	--min-MQ 60 \
	--remove-tag OQ \
	--input-fmt-option 'filter=ncigar==1 && ![XA] && ![SA]' \
	"${input_bam}"

samtools index "${output_bam}"

if [ ${pipe_echo} -eq 1 ]; then echo "${output_bam}"; fi

