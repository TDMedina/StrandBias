#!/bin/bash
#SBATCH --partition=normal,highmem,MSC
#SBATCH --job-name=flagstat

#SBATCH --time=00-00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

#SBATCH --error=%x.%j.stderr
#SBATCH --output=%x.%j.stdout

#SBATCH --mail-user=tyler.medina@universityofgalway.ie
#SBATCH --mail-type=BEGIN,END,FAIL

# Start Script

function help {
	cat <<- EOF

	00.flagstat.sh [-o <output.flagstat>] -i <input.bam>

	-i <input.bam>	Input BAM file.
	-o <output.flagstat>	Output flagstat file. [Default=<input.bam>.flagstat]

	EOF
}

while getopts ":i:o:h" arg; do
	case "${arg}" in
		h)
			help
			exit 0
			;;
		i)
			input_bam="${OPTARG}"
			;;
		o)
			output_file="${OPTARG}"
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

if [ -z ${output_file+x} ]; then output_file="${input_bam}.flagstat"; fi

samtools flagstat -O tsv "${input_bam}" > "${output_file}"
