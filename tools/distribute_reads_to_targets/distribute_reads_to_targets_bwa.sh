#!/usr/bin/env bash

# sanity check
printf "Conda env: $CONDA_DEFAULT_ENV\n"
printf "Python version: $(python --version |  awk '{print $2}')\n"
#printf "Python version: $(conda list | egrep python3 | awk '{print $2}')\n"
printf "Biopython version: $(conda list | egrep biopython | awk '{print $2}')\n"
printf "Samtools version: $(conda list | egrep samtools | awk '{print $2}')\n"
printf "Unzip version: $(unzip -v | head -n1 | awk '{print $2}')\n"
printf "Zip version: $(zip -v | head -n1 | awk '{print $2}')\n"
printf "Bash version: ${BASH_VERSION}\n\n"

# The runDistributeToTargets function calls the python script, which,
# after a BWA search against the target sequences, sorts the hits
# from the BAM-file to the successful hits.
# Then copies the output to one final ZIP file and removes the temporary files.

# Usage: sh distribute_reads_to_targets_bwa.sh -b Galaxy82-[Merged_BAM].bam -r Raw_reads_test.zip -o distributed_reads.zip
#

runDistributeToTargets() {
    strScriptDir=$(dirname "$(readlink -f "$0")")
    base_location=$(echo ${outputzip} | egrep -o '^.*files')
    strDirectory=$(mktemp -d ${base_location}/XXXXXX)
    mkdir -p "${strDirectory}_temp"
    mkdir -p "${strDirectory}_temp/merged_reads"
    mkdir -p "${strDirectory}_temp/outputfile/"
    printf "Base_location: ${base_location}\n" >&2
    printf "strDirectory location: ${strDirectory}\n" >&2
    unzip ${readfolder} -d ${strDirectory}_temp/packed_reads
    unzip ${strDirectory}_temp/packed_reads"/*.zip" -d ${strDirectory}_temp/raw_reads

# Merge de FASTQ files van elke pair in twee losse fastq files
    cat ${strDirectory}_temp/raw_reads/*R1*.fastq > ${strDirectory}_temp/merged_reads/merged_R1.fastq
    cat ${strDirectory}_temp/raw_reads/*R2*.fastq > ${strDirectory}_temp/merged_reads/merged_R2.fastq

# Gebruikt de de unpacked fastq files als input en outputs (losse fasta files in de -o flag)
    python3 $strScriptDir"/distribute_reads_to_targets_bwa.py" -b ${bamfile} \
                                               -r ${strDirectory}_temp/merged_reads \
                                               -o ${strDirectory}_temp/fastafiles/
    printf "Python output: $(ls ${strDirectory}_temp/fastafiles/)\n" >&2
    printf "Distributing script finished successfully, attempting to create Zip...\n" >&2

# Zip command die de losse fasta files uit het python script in een zip bestand doet in een tijdelijke folder
    zip -rj ${strDirectory}_temp/tempzip.zip ${strDirectory}_temp/fastafiles/gene*

# Kopieer de gezipte fasta files naar de uiteindelijke outputlocatie gegenereert door Galaxy
    cp ${strDirectory}_temp/tempzip.zip ${outputzip}
    printf "_temp folder contents: $(ls ${strDirectory}_temp)\n" >&2
#    cat ${strDirectory}_temp/tempzip.zip > ${outputzip}
#    mv ${strDirectory}_temp/tempzip.zip ${outputzip}
    rm -rf ${strDirectory}_temp
    printf "Outputlocation: "${outputzip}"\n" >&2
    printf "Output.dat folder contents: "$(ls ${outputzip})"\n" >&2
    printf "Shell script finished successfully\n" >&2

}

# The main function.
main() {
    runDistributeToTargets
}

# The getopts function.
while getopts ":b:r:o:vh" opt; do
    case ${opt} in
        b)
            bamfile=${OPTARG}
            ;;
        r)
            readfolder=${OPTARG}
            ;;
        o)
            outputzip=${OPTARG}
            ;;
        v)
            echo ""
            echo "distribute_reads_to_targets_bwa.sh [1.6.3]"
            echo ""

            exit
            ;;
        h)
            echo ""
            echo "Usage: distribute_reads_to_targets_bwa.sh [-h] [-v] "
            echo "                 [-b INPUT_BAM] [-r RAW_READS_ZIP]  "
            echo "                 [-o OUTPUT_ZIP]                    "
            echo ""
            echo "After a BWA search of the raw reads against "
            echo "the target sequences, the reads need to be sorted "
            echo "according to the successful hits. "
            echo "This script takes the BWA output (BAM format) "
            echo " and the raw read files"
            echo "and distributes the reads into FASTA files "
            echo "ready for assembly."
            echo ""

            exit
            ;;
        \?)
            echo ""
            echo "You've entered an invalid option: -${OPTARG}."
            echo "Please use the -h option for correct formatting information."
            echo ""

            exit
            ;;
        :)
            echo ""
            echo "You've entered an invalid option: -${OPTARG}."
            echo "Please use the -h option for correct formatting information."
            echo ""

            exit
            ;;
    esac
done

main

# Additional information:
# =======================
#
# If there are multiple results (for example, one for each read direction),
# concatenate them prior to sorting.