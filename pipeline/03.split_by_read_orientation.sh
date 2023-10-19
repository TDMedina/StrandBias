#!/bin/bash
#SBATCH --partition=normal,highmem,MSC
#SBATCH --job-name=split_by_read_orientation


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
		
		03.split_by_read_orientation.sh [-r <reference.fasta> -o <output.bam>] -d {F1,F2,R1,R2} -i <input.bam>
		
		-i <input.bam>	Input BAM file.
		-o <output.bam>	Output BAM file. [Default=<input>.<orienation>.bam]
		-d {F1,F2,R1,R2}	Read orientation to extract.
		-r <reference.fasta>	Reference FASTA file. [Default=./reference.fa]
		-p	If set, script will echo the output name for pipeline purposes.		
	EOF
}

reference="./reference.fa"
pipe_echo=0

while getopts ":r:i:o:d:ph" arg; do
	case "${arg}" in
		r)
			reference="${OPTARG}"
			;;
		i)
			input_bam="${OPTARG}"
			;;
		o)
			output_bam="${OPTARG}"
			;;
		d)
			orientation="${OPTARG}"
			case "${OPTARG}" in
				F1)
					flag=96
					;;
				F2)
					flag=160
					;;
				R1)
					flag=80
					;;
				R2)
					flag=144
					;;
				*)
					echo "Invalid value for option -d; Must be one of {F1, F2, R1, R2}"
					help
					exit 1
					;;
			esac
			;;
		p)
			pipe_echo=1
			;;
		h)
			help
			exit 0
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

if [ -z ${input_bam+x} ]; then echo "Missing required option: -i"; help; exit 1; fi
if [ -z ${orientation+x} ]; then echo "Missing required option: -d"; help; exit 1; fi

if [ -z ${output_bam+x} ]; then output_bam="${input_bam%%.bam}.${orientation}.bam"; fi


# echo "input=${input_bam} output=${output_bam} orientation=${orientation} flag=${flag}"

samtools view \
	--reference "${reference}" \
	--output "${output_bam}" \
	--bam \
	--require-flags "${flag}" \
	"${input_bam}"

if [ ${pipe_echo} -eq 1 ]; then echo "${output_bam}"; fi
