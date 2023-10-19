#!/bin/bash
#SBATCH --partition=normal,highmem,MSC
#SBATCH --job-name=pileup


#SBATCH --time=00-00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

#SBATCH --error=%x.%j.stderr
#SBATCH --output=%x.%j.stdout

#SBATCH --mail-user=tyler.medina@universityofgalway.ie
#SBATCH --mail-type=BEGIN,END,FAIL

# Start Script

function help {
	cat <<- EOF

	04.pileup.sh [-r <reference.fa> -o <output.pileup>] -i <input.bam>

	-i <input.bam>	Input BAM file.
	-o <output.pileup>	Output pileup file. [Default=<input>.pileup]
	-r <reference.fasta>	Reference FASTA file. [Default=./reference.fa]
	-q <int>	Minimum BQ score required [Default=37]

	EOF
}

reference="./reference.fa"
bq=37

while getopts ":r:i:o:q:h" arg; do
	case "${arg}" in
		h)
			help
			exit 0
			;;
		i)
			# input_bams="${OPTARG}"
			mapfile -d "," -t input_bams < <(echo -n "${OPTARG}")
			;;
		o)
			output_file="${OPTARG}"
			;;
		r)
			reference="${OPTARG}"
			;;
		q)
			bq="${OPTARG}"
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
if [ -z ${input_bams+x} ]; then echo "Missing required argument: -i"; help; exit 1; fi

if [ -z ${output_file+x} ]; then output_file="${input_bams[0]%%.bam}.pileup"; fi

samtools mpileup \
	--fasta-ref "${reference}" \
	--no-BAQ \
	--min-BQ "${bq}" \
	--output-BP \
	--reverse-del \
	"${input_bams[@]}" > "${output_file}" 2> "${output_file%.*}.stderr"
#	"${input_bam}" > "${output_file}" 2> "${output_file%.*}.stderr"

awk '$5 !~ /^((\^.)?([\.,*#])(\$)?)+$/ {print}' "${output_file}" > "${output_file%.pileup}.no_match_positions.pileup"
