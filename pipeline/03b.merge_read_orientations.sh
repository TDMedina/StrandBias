#!/bin/bash
#SBATCH --partition=normal,highmem,MSC
#SBATCH --job-name=merge_read_orientations


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
	
	03b.merge_read_orientations.sh [-r <reference.fa> -o <output.bam> -p] -a <input1.bam> -b <input1.bam> -d <orientation_label>

	-a <input1.bam>	First input BAM file.
	-b <input2.bam>	Second input BAM file.
	-d <orientation_label>	Label to use for default output name. Ignored if -o is used.
	-o <output.bam>	Output BAM file. [Default=<input1>.<orientation>.bam]
	-r <reference.fasta>	Reference FASTA file. [Default=./reference.fa]
	-p	If set, script will echo the output name for pipeline purposes.

	EOF
}

reference="./reference.fa"
pipe_echo=0

while getopts ":r:a:b:d:o:ph" arg; do
	case "${arg}" in
		h)
			help
			exit 0
			;;
		a)
			input_bam_1="${OPTARG}"
			;;
		b)
			input_bam_2="${OPTARG}"
			;;
		d)
			label="${OPTARG}"
			;;
		o)
			output_bam="${OPTARG}"
			;;
		r)
			reference="${OPTARG}"
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
if [ -z ${input_bam_1+x} ]; then echo "Missing required argument: -a"; help; exit 1; fi 
if [ -z ${input_bam_2+x} ]; then echo "Missing required argument: -b"; help; exit 1; fi 

if [ -z ${output_bam+x} ]; then output_bam="${input_bam%%.bam}.filtered.bam"; fi 
if [ -z ${output_bam+x} ]; then if [ -z ${label+x} ]; then echo "Missing required argument for default output naming: -d"; help; exit 1; fi; output_bam="${input_bam_1%.*.bam}.${label}.bam"; fi

samtools merge \
	-o "${output_bam}" \
	"${input_bam_1}" \
	"${input_bam_2}"

if [ ${pipe_echo} -eq 1 ]; then echo "${output_bam}"; fi
