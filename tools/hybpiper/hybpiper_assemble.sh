#!/bin/bash
# #!/usr/bin/env bash

# sanity check
#printf "Conda env: $CONDA_DEFAULT_ENV\n" >&1
#printf "Biopython version: $(conda list | egrep biopython | awk '{print $2}')\n" >&1
#printf "Samtools version: $(conda list | egrep samtools | awk '{print $2}')\n" >&1
#printf "Unzip version: $(unzip -v | head -n1 | awk '{print $2}')\n" >&1
#printf "Zip version: $(zip -v | head -n1 | awk '{print $2}')\n" >&1
#printf "Bash version: ${BASH_VERSION}\n\n" >&1

# The runDistributeToTargets function calls the python script, which,
# after a BWA search against the target sequences, sorts the hits
# from the BAM-file to the successful hits.
# Then copies the output to one final ZIP file and removes the temporary files.

# Usage: sh distribute_reads_to_targets_bwa.sh -b Galaxy82-[Merged_BAM].bam -r Raw_reads_test.zip -o distributed_reads.zip
#

# version 1.0.0 (nonfunctional)

runHybpiperAssemble() {
    (
    strScriptDir=$(dirname "$(readlink -f "$0")")
#    strScriptDir="/home/eremus007/Desktop/Hybpiper/scripts"
    base_location=$(echo ${outputzip} | egrep -o '^.*files')
#    base_location="/home/eremus007/Desktop/Hybpiper/final_output/"
    strDirectory=$(mktemp -d ${base_location}/XXXXXX)

    mkdir -p "${strDirectory}_temp"
    mkdir -p "${strDirectory}_temp/outputfile/"
    outputDirectory=$(dirname "${readfiles}")


    #Unzip readfile zip and copy to proper location
    unzip ${readfiles} -d ${strDirectory}_temp/raw_reads
    rawReadsFile="${strDirectory}_temp/raw_reads"

    # Generate Hybpiper commands
     python3 ${strScriptDir}/generate_hybpiper_commands.py \
        -r ${rawReadsFile} \
        -o "cmdfile.txt" \
        -t ${targetfile} \
        -f ${target_format} \
        -e ${search_engine} \
        -i ${intronerate_bool} \
        -m ${heatmap_bool}

    # Execute generated Hybpiper commands
    while read cmd_to_execute
      do
        eval "${cmd_to_execute}"
      done < cmdfile.txt

    # move everything to proper location
    cp -r ${outputDirectory}/* ${strDirectory}_temp/outputfile/

    # Zip the correct output files in the temporary directory
    zip -r ${strDirectory}_temp/tempzip.zip ${strDirectory}_temp/outputfile/*
    # Copy the zip file to the galaxy data.dat
    cp ${strDirectory}_temp/tempzip.zip ${outputzip}

    # Delete remaining temporary files
    rm -rf ${strDirectory}_temp

    # Delete Generated commands txt file
    rm ${strScriptDir}/cmdfile.txt

    ) #> /dev/null 2>&1 #This will make sure nothing is written to stderr (hopefully)
}

# The main function.
main() {
  runHybpiperAssemble
}

# The getopts function.
while getopts ":r:o:t:f:e:i:m:vh" opt; do
    case ${opt} in
        r)
            readfiles=${OPTARG}
            ;;
        o)
            outputzip=${OPTARG}
            ;;
        t)
            targetfile=${OPTARG}
            ;;
        f)
            target_format=${OPTARG}
            ;;
        e)
            search_engine=${OPTARG}
            ;;
        i)
            intronerate_bool=${OPTARG}
            ;;
        m)
            heatmap_bool=${OPTARG}
            ;;
        v)
            echo ""
            echo "hybpiper_assemble.sh [x.x.x]"
            echo ""

            exit
            ;;
        h)
            echo ""
            echo "Usage: distribute_reads_to_targets_bwa.sh [-h] [-v] "
            echo "                 [-r RAW_READS_ZIP]                 "
            echo "                 [-o OUTPUT_ZIP]                    "
            echo "                 [-t SEQUENCE_TYPE]                   "
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
