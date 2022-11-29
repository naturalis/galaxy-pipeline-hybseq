#!/usr/bin/env bash

# sanity check
#printf "Conda env: $CONDA_DEFAULT_ENV\n"
#printf "Python version: $(python --version |  awk '{print $2}')\n"
#printf "Biopython version: $(conda list | egrep biopython | awk '{print $2}')\n"
#printf "Samtools version: $(conda list | egrep samtools | awk '{print $2}')\n"
#printf "Unzip version: $(unzip -v | head -n1 | awk '{print $2}')\n"
#printf "Bash version: ${BASH_VERSION}\n\n"

# The runDistributeToTargets function calls the python script, which,
# after a BWA search against the target sequences, sorts the hits
# from the BAM-file to the successful hits.
# Then copies the output to one final ZIP file and removes the temporary files.
#

runDistributeToTargets() {
    strScriptDir=$(dirname "$(readlink -f "$0")")
    strDirectory=$(mktemp -d)
    mkdir -p "${strDirectory}_temp"
    unzip ${strScriptDir}/${readfolder} -d ${strScriptDir}/packed_reads
    unzip ${strScriptDir}/packed_reads"/*.zip" -d ${strScriptDir}/raw_reads
    rm -rf packed_reads
    python3 $strScriptDir"/distribute_reads_to_targets_bwa.py" -b ${bamfile} \
                                               -r $strScriptDir/raw_reads
    zip -r ${strScriptDir}/${outputzip} gene*
    rm -rf gene*
    rm -rf raw_reads

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
            echo "distribute_reads_to_targets_bwa.sh [1.4.0]"
            echo ""

            exit
            ;;
        h)
            echo ""
            echo "Usage: distribute_reads_to_targets_bwa.sh [-h] [-v] "
            echo "                 [-b INPUT_BAM] [-r RAW_READS_FOLDER]    "
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
