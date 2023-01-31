#!/usr/bin/env python

'''
Version: 1.0.1
Author: Jeremy van Veen
'''
import os, argparse


def get_file_names(directory):
    """Iterates through a specified directory and returns a list with
    the file names.

        Parameters
        ----------
        directory : str
            The file containing the files to iterate through

        Returns
        -------
        filenames : list
            a list of strings which are the names of the files in the
            specified directory.
        """
    # iterate over files in
    # that directory
    filenames = []
    for filename in os.listdir(directory):
        file = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(file):
    # split off extentions
            file_list = file.split('.')
    # remove directory path/extentions and keep file name
            filenames.append(file_list[0][1+len(directory):])
    return filenames


def get_sample_names(filenames, unique=False):
    samplenames_list, pair_dict = [], {}
    for filename in filenames:
        part_index = 0
        filename_parts = filename.split("_")
        if "R1" in filename_parts:
            part_index = filename_parts.index("R1")
        elif "R2" in filename_parts:
            part_index = filename_parts.index("R2")
        samplename = filename_parts[:part_index][0]
        samplenames_list.append(samplename)
        if unique:
            samplenames_list = set(samplenames_list)
            samplenames_list = list(samplenames_list)
            samplenames_list.sort()
    return samplenames_list

def write_commands_to_file(to_write_list, outputlocation="cmdfile.txt", overwrite=False):
    file = outputlocation
    if overwrite:
        resetfile = open(file, "w")
        resetfile.close()
    with open(file, "a") as outfile:
        for element in to_write_list:
            outfile.write(element+"\n")


def construct_assemble_commands(readfile, samplenames, target_format, targetfile, search_engine, intronerate_bool):
    assemble_cmds = []
    for sample in samplenames:
        assemble_cmd = "hybpiper assemble -r %s/%s.fastq" % (str(readfile), sample)
        if target_format == "AA" or target_format == "aa":
            assemble_cmd= str(assemble_cmd) + " -t_aa %s" % targetfile
        elif target_format == "DNA" or target_format == "dna":
            assemble_cmd = str(assemble_cmd) + " -t_dna %s" % targetfile
        assemble_cmd = str(assemble_cmd) + " --prefix %s" % sample

        if search_engine == "diamond":
            assemble_cmd = str(assemble_cmd) + " --diamond"
        elif search_engine == "bwa":
            assemble_cmd = str(assemble_cmd) + " --bwa"

        if intronerate_bool == "y":
            assemble_cmd = str(assemble_cmd) + " --run_intronerate"
        # assemble_cmd = str(assemble_cmd) + " --hybpiper_dir %s" % str(output_location)

        assemble_cmds.append(assemble_cmd)
    return assemble_cmds


def construct_stats_command():
    stat_cmd = "hybpiper stats -t_dna test_targets.fasta gene namelist.txt"
    return stat_cmd

def construct_heatmap_command():
    heatmap_cmd = "hybpiper recovery_heatmap seq_lengths.tsv"
    return heatmap_cmd


def construct_retrieve_commands():
    retrieve_cmds = []
    retrieve_cmd_1 = "hybpiper retrieve_sequences -t_dna test_targets.fasta dna --sample_names namelist.txt --fasta_dir 01_dna_seqs"
    retrieve_cmd_2 = "hybpiper retrieve_sequences -t_dna test_targets.fasta aa --sample_names namelist.txt --fasta_dir 02_aa_seqs"
    retrieve_cmds.append(retrieve_cmd_1)
    retrieve_cmds.append(retrieve_cmd_2)

    return retrieve_cmds

def construct_paralog_command():
    paralog_cmd = "hybpiper paralog_retriever namelist.txt -t_dna test_targets.fasta"
    return paralog_cmd


def parseArgvs():
    parser = argparse.ArgumentParser(prog="generate_namelist.py",
                                     description="Script to create a file "
                                                 "containing the name of "
                                                 "every sample in the readfile ")
    parser.add_argument("-v", "--version", action="version",
                        version="generate_namelist.py 1.0.1")
    parser.add_argument("-r", "--readfile", action="store", dest="readfile",
                        help="The location of the input readfile(s)",
                        required=True)
    parser.add_argument("-o", "--output", action="store", dest="output_path",
                        help="The path location where the output file with the new commands should go",
                        required=False)
    parser.add_argument("-t", "--targets", action="store", dest="targetfile",
                        help="The location of the input targetfile",
                        required=True)
    parser.add_argument("-f", "--target_format", action="store", dest="target_format",
                        help="The type of sequences the targets are comprised off, either DNA or Amino acids.",
                        required=True)
    parser.add_argument("-e", "--engine", action="store", dest="search_engine",
                        help="Which method Hybpiper uses to map the sequences to the targets.", required=True)
    parser.add_argument("-i", "--intronerate", action="store", dest="intronerate_bool",
                        help="Whether Hybpiper should extract the introns using intronerate", required=False)
    parser.add_argument("-m", "--heatmap_bool", action="store", dest="heatmap_bool",
                        help="Whether hybpiper should generate a heatmap of the gene recovery", required=False)
    parser.add_argument("-g", "--output_g", action="store", dest="output_galaxy",
                        help="The path location where the output should go for Galaxy",
                        required=False)
    argvs = parser.parse_args()
    return argvs


def main():
    argvs = parseArgvs()
    readfile = argvs.readfile
    targetfile = argvs.targetfile
    target_format = argvs.target_format
    search_method = argvs.search_engine
    intronerate_bool = argvs.intronerate_bool
    galaxy_output_location = argvs.output_galaxy
    heatmap_bool = argvs.heatmap_bool
    output_location = argvs.output_path

    filenameslist = get_file_names(readfile)
    samplenameslist = get_sample_names(filenameslist, True)
    # filenameslist = get_sample_names(filenames, True)
    ## Write samplenames to namelist.txt
    # if output_location == None:
    #     write_Samplenames_to_file(samplenameslist)
    # else:
    #     write_Samplenames_to_file(samplenameslist, output_location)
    cmds = construct_assemble_commands(readfile, filenameslist, target_format, targetfile, search_method, intronerate_bool)

    stats_cmd = construct_stats_command()
    heatmap_cmd = construct_heatmap_command()
    retrieve_cmds = construct_retrieve_commands()
    paralog_cmd = construct_paralog_command()


    # write commands to .txt file
    write_commands_to_file(cmds, output_location, overwrite=True)
    #write_commands_to_file([stats_cmd])
    #write_commands_to_file([heatmap_cmd])
    #write_commands_to_file(retrieve_cmds)
    #write_commands_to_file([paralog_cmd])

if __name__ == '__main__':
    main()
