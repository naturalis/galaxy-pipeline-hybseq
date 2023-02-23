#!/bin/bash

# The runHybpiperAssemble function calls the python script, which,
# generates the nessesary hybpiper commands based on the inputs and writes
# them to a .txt file.
# This text file is then opened in a while loop and executes every command
# one by one.
#
# Usage:
# sh hybpiper_assemble.sh -r <readfiles.zip> -o <NAME_FOR_OUTPUT.zip> -t <test_targets.fasta> -f <dna/aa> -e <bwa/diamond/default> -i <y/n> -m <y/n>
#
#

runHybpiperExtractLoci() {
    (
    strScriptDir=$(dirname "$(readlink -f "$0")")
#    strScriptDir="/home/eremus007/Desktop/Hybpiper/scripts"
    # base_location=$(echo ${outputfolder} | egrep -o '^.*files')
    base_location=$(dirname "${outputfolder}")
#    base_location="/home/eremus007/Desktop/Hybpiper"
    strDirectory=$(mktemp -d "${base_location}"/XXXXXX_temp)
    workingDir="${strDirectory}/working_dir"

    # Create Temporary Directories
    mkdir -p "${workingDir}"
    mkdir -p "${strDirectory}/outputfile/"
    mkdir -p "${workingDir}/hybpiper_output"
    mkdir -p "${workingDir}/extracted_loci"

    # Change the working directory to be a temporary folder
    cd "${workingDir}" || exit

    # Unzip hybpiper output zip
    unzip "${inputfolder}" -d "${workingDir}"/hybpiper_output

    # Extract the loci sequences using python script
     python3 "${strScriptDir}"/extract_loci_sequences.py \
        -f "${workingDir}"/hybpiper_output \
        -o "${workingDir}"/extracted_loci \
        -t "${file_type}"

    # Zip the correct output files in the temporary directory
    7z a "${strDirectory}"/tempzip.zip "${workingDir}"/extracted_loci/*

    # Copy the zip file to the galaxy data.dat
    cp "${strDirectory}"/tempzip.zip "${outputfolder}"

    # Delete remaining temporary files
    rm -rf "${strDirectory}"

    ) > /dev/null 2>&1 #This will make sure nothing is written to stderr
}

# The main function.
main() {
  runHybpiperExtractLoci
}

# The getopts function.
while getopts ":f:o:t:vh" opt; do
    case ${opt} in
        f)
            inputfolder=${OPTARG}
            ;;
        o)
            outputfolder=${OPTARG}
            ;;
        t)
            file_type=${OPTARG}
            ;;
        v)
            echo ""
            echo "extract_loci_sequences.sh [0.1.1]"
            echo ""

            exit
            ;;
        h)
            echo ""
            echo "Usage: sh extract_loci_sequences.sh [-h] [-v]       "
            echo "                 [-f HYBPIPER_OUTPUT_ZIP.zip        "
            echo "                 [-o OUTPUT_ZIP.zip]                "
            echo "                 [-t FILE_TYPE]                     "
            echo ""
            echo "HybPiper was designed for targeted sequence capture,"
            echo "in which DNA sequencing libraries are enriched for gene regions of interest, especially for phylogenetics. "
            echo "HybPiper is a suite of Python scripts/modules that wrap and connect bioinformatics tools "
            echo "to extract target sequences from high-throughput DNA sequencing reads "
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