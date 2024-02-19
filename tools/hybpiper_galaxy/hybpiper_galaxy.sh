#!/bin/bash
set -e
# The runHybpiperAssemble function calls the python script, which,
# generates the necessary hybpiper commands based on the inputs and writes
# them to a .txt file.
# This text file is then opened in a while loop and executes every command
# one by one.
#
# Usage:
# sh hybpiper_galaxy.sh -r <readfiles.zip> -o <NAME_FOR_OUTPUT.zip> -t <test_targets.fasta> -f <dna/aa> -e <bwa/diamond/default> -i <y/n> -m <y/n>
#
#

runHybpiperAssemble() {
    (
    strScriptDir=$(dirname "$(readlink -f "$0")")
    base_location=$(dirname ${outputzip}) 
    #base_location="/home/eremus007/Desktop/Hybpiper" #Residual for testing outside of Galaxy, might still come in handy
    strDirectory=$(mktemp -d ${base_location}/XXXXXX_temp)
    workingDir="${strDirectory}/working_dir"

    mkdir -p ${workingDir}
    #mkdir -p "${strDirectory}/outputfile/"
    mkdir -p "${workingDir}/raw_reads"

    #cd ${workingDir}

    #outputDirectory=$(dirname "${readfiles}")


    #Unzip readfile zip and copy to proper location
    cp ${targetfile} ${workingDir}/targetfile.fasta
    # unzip ${readfiles} -d ${workingDir}/raw_reads
    rawReadsFile="${workingDir}/raw_reads"

    # Generate Hybpiper commands
    no_intronerate_arg=$([ "$run_intronerate" = "n" ] && echo "--no_intronerate" || echo "")
    heatmap_arg=$([ "$heatmap" = "y" ] && echo "--heatmap" || echo "")
    echo $no_intronerate_arg
    echo $heatmap_arg
    python3 ${strScriptDir}/generate_hybpiper_commands.py \
        --read_zip "${readfiles}" \
        --read_dir "${rawReadsFile}" \
        --output_dir "${workingDir}" \
        --targets "${workingDir}/targetfile.fasta" \
        --target_format "${target_format}" \
        --engine "${search_engine}" \
        ${no_intronerate_arg} \
        ${heatmap_arg} \
        --write_namelist \
        --timeout "${timeout}"

    # Execute generated Hybpiper commands
    while read cmd_to_execute
      do
        eval "${cmd_to_execute}"
      done < "${workingDir}/cmdfile.txt"

    # move everything to proper location but before that remove the copies of the input files so they don't get copied with
    #rm -r ${rawReadsFile}
    #rm ${workingDir}/targetfile.fasta
    #cp -rf ${workingDir}/* ${strDirectory}/outputfile/

    # Zip the correct output files in the temporary directory
    7z a ${strDirectory}/tempzip.zip ${workingDir}/assemble/*
    # Copy the zip file to the galaxy data.dat
    cp ${strDirectory}/tempzip.zip ${outputzip}

    # Delete remaining temporary files
    #rm -rf ${strDirectory}

    ) 2> /home/galaxy/error.txt
}

# The main function.
main() {
  runHybpiperAssemble
}

# The getopts function.
while getopts ":r:o:t:f:e:i:m:x:vh" opt; do
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
            run_intronerate=${OPTARG}
            ;;
        m)
            heatmap=${OPTARG}
            ;;
        x)
            timeout=${OPTARG}
            ;;
        v)
            echo ""
            echo "hybpiper_galaxy.sh [1.2.0]"
            echo ""

            exit
            ;;
        h)
            echo ""
            echo "Usage: hybpiper_galaxy.sh [-h] [-v] "
            echo "                 [-r RAW_READS_ZIP]                 "
            echo "                 [-o OUTPUT_ZIP]                    "
            echo "                 [-t TARGET_FILE_FASTA]             "
            echo "                 [-f TARGET_FORMAT]                 "
            echo "                 [-e MAPPING METHOD]                "
            echo "                 [-i RUN_INTRONERATE]               "
            echo "                 [-m HEATMAP]                       "
            echo "                 [-x TIME_OUT_INT]                  "
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
