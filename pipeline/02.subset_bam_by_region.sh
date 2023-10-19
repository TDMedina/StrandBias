#!/bin/bash
#SBATCH --partition=normal,highmem,MSC
#SBATCH --job-name=subset_bam_by_region


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

	02.subset_bam_by_region.sh [-p -r <reference.fa>] [-o <output.bam> | -s <region_name>] -i <input.bam> -b <regions.bed>

	-i <input.bam>	Input BAM file.
	-b <regions.bed>	BED file of desired regions.
	-s <region_name>	Label to use for default output name. Ignored if -o is used.
	-o <output.bam>	Output BAM file. [Default=<input>.<region_name>.bam]
	-r <reference.fasta>	Reference FASTA file. [Default=./reference.fa]
	-p	If set, script will echo the output name for pipeline purposes.	

	EOF
}

reference="./reference.fa"
pipe_echo=0

while getopts ":i:b:s:o:r:ph" arg; do
	case "${arg}" in
		h)
			help
			exit 0
			;;
		i)
			input_bam="${OPTARG}"
			;;
		b)
			bed="${OPTARG}"
			;;
		s)
			region_name="${OPTARG}"
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
if [ -z ${input_bam+x} ]; then echo "Missing required argument: -i"; help; exit 1; fi
if [ -z ${bed+x} ]; then echo "Missing required argument: -b"; help; exit 1; fi

if [ -z ${output_bam+x} ]; then if [ -z ${region_name+x} ]; then echo "Missing required argument for default output naming: -s"; help; exit 1; fi; output_bam="${input_bam%%.bam}.${region_name}.bam"; fi


samtools view \
	--reference "${reference}" \
	--region-file "${bed}" \
	--output "${output_bam}" \
	--bam \
	"${input_bam}"

samtools index "${output_bam}"

if [ ${pipe_echo} -eq 1 ]; then echo "${output_bam}"; fi
