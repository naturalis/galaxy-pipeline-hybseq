#!/bin/bash

# The runHybpiperAnalysis function calls the python script, which,
# concatenates the input .fa files into one .fa file per locus.
# Then, runs MUSCLE seperately for each concatenated .fa file
# to generate alignment files in phylip (.phy) format.
# Those alignments are then concatenated into a supermatrix
# which itself is also in .phy format alongside a .txt file with the 
# partition information of each loci in the supermatrix.
#
#
# Usage:
# sh hybpiper_analysis.sh -i <EXTRACTED_LOCI_SEQUENCES.zip> -o <NAME_FOR_OUTPUT.zip>
#
#

runHybpiperAnalysis() {
    (
    strScriptDir=$(dirname "$(readlink -f "$0")")
    base_location=$(dirname "${outputfolder}")

    strDirectory=$(mktemp -d "${base_location}"/XXXXXX_temp)
    workingDir="${strDirectory}/working_dir"

    # Create Temporary Directories
    mkdir -p "${workingDir}"
    mkdir -p "${strDirectory}/outputfile/"
    mkdir -p "${workingDir}/analysis_output"

    # Change the working directory to be a temporary folder
    cd "${workingDir}" || exit

    # Unzip hybpiper output zip
    unzip "${inputfolder}" -d "${workingDir}"/extracted_loci

    # Extract the loci sequences using python script
     python3 "${strScriptDir}"/hybpiper_analysis.py \
        -i "${workingDir}"/extracted_loci  \
        -o "${workingDir}"/analysis_output

    # Zip the correct output files in the temporary directory
    7z a "${strDirectory}"/tempzip.zip "${workingDir}"/analysis_output/*

    # Copy the zip file to the galaxy data.dat
    cp "${strDirectory}"/tempzip.zip "${outputfolder}"

    # Delete remaining temporary files
    rm -rf "${strDirectory}"

    ) #> /dev/null 2>&1 #This will make sure nothing is written to stderr
}

# The main function.
main() {
  runHybpiperAnalysis
}

# The getopts function.
while getopts ":i:o:vh" opt; do
    case ${opt} in
        i)
            inputfolder=${OPTARG}
            ;;
        o)
            outputfolder=${OPTARG}
            ;;
        v)
            echo ""
            echo "hybpiper_analysis.sh [1.0.2]"
            echo ""

            exit
            ;;
        h)
            echo ""
            echo "Usage: sh hybpiper_analysis.sh [-h] [-v]       "
            echo "                 [-i EXTRACTER_LOCI_ZIP.zip    "
            echo "                 [-o OUTPUT_ZIP.zip]           "
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