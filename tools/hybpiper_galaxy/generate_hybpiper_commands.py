#!/usr/bin/env python

'''
This script takes several inputs using the argparse command, with which it
generates the commands needed to run Hybpiper.

It works in conjunction with- and is called by the hybpiper_galaxy.sh shell
script to create a file with the complete commands.

The shell script would then go over each command in that text file one by one
and execute them to run the Hybpiper pipeline as intended.
'''
import os, argparse


def get_file_names(directory, remove_file_extentions=True):
    """Iterates through a specified directory and returns a list with
    the file names (without the extention)
    Useful for getting the names of the sample.

        Parameters
        ----------
        directory : str
            The file containing the files to iterate through
        remove_file_extentions : bool
            A boolean indicating whether the function should remove
            file extentions or keep them.

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
            # Make sure all slash charactersare the same
            file = file.replace('\\', '/')
            # split off extentions
            file_list = file.split('/')
            extracted_filename = file_list[-1:][0]
            # remove directory path/extentions and keep file name if
            # boolean is not true
            if remove_file_extentions:
                cut_filename = extracted_filename.split('.')
                filenames.append(cut_filename[0])
            else:
                filenames.append(extracted_filename)
    return filenames


def get_sample_names(filenames, unique=False):
    """Iterates through a list of filenames, splits the filenames on every
    underscore they contain, and appends only the part before the first
    underscore to the list it returns.
    The output is a list containing just the sample name, and none of the other
    extensions like _R1 or _test

            Parameters
            ----------
            filenames : list
                The list with all the full filenames
            unique : bool (default=False)
                A boolean indicating whether the output list should only
                have unique values.
                If true, the function removes duplicates by casting the list
                as a set and then back to a list.

            Returns
            -------
            samplenames_list : list
                a list of strings which are just the sample names
                instead of the full file name.
        """
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


def write_commands_to_file(to_write_list, outputlocation=".", outputfilename = 'cmdfile.txt',
                           overwrite=False):
    """Iterates through a list, and writes each element to a new line in the
    specified output file.

                Parameters
                ----------
                to_write_list : list
                    A list with all the elements the function needs to write
                    to the output file.
                outputlocation : str
                    A string with the path to where the outputfile should be
                    by default the output directory is the current
                    working directory
                outputfilename : str
                    A string with the path to where the outputfile should be
                    by default the output file is a txt file called
                    'cmdfile.txt' in the same directory as the script.
                overwrite : bool
                    A boolean indicating whether the the outputfile should be
                    reset. If false, output is apended to the output file, if
                    True, the output file is emptied before being written to.
    """
    file = str(outputlocation) + '/%s' % outputfilename
    if overwrite:
        resetfile = open(file, "w")
        resetfile.close()
    with open(file, "a") as outfile:
        for element in to_write_list:
            outfile.write(element + "\n")


def construct_assemble_commands(readfile, filenames, samplenames,
                                target_format, targetfile,
                                search_engine, intronerate_bool):
    """The function constructs the hybpiper assemble command
    It does this by appending a template string with the proper flags according
    to the state of the arguments. Using this method, it assembles one
    command for each sample for every read file. Then returns a list
    of all the generated commands.

                Parameters
                ----------
                readfile : str
                    A string resembling the path to the folder with the read
                    fastq file.
                filenames : list
                    A list containing the names of all the files.
                samplenames : list
                    A list containing the names of all the samples (different
                    from filenames).
                target_format : str
                    The list with all the full filenames
                targetfile : str
                    A string resembling the path to the targets fasta file.
                search_engine : str
                    A string coresponding to the type of method hybpiper
                    should use for mapping reads.
                intronerate_bool : str
                    Functions as a boolean, with values 'y' and 'n'
                    to represent true and false respectively. This 'boolean'
                    checks to see to add the flag to run intronerate
                    to the hybpiper assemble command.

                Returns
                -------
                assemble_cmds : list
                    a list of the generated hybpiper assemble commands in the
                    form of strings.
            """
    assemble_cmds = []
    for file in filenames:
        assemble_cmd = "hybpiper assemble -r %s/%s" % (
        str(readfile), file)
        if target_format == "AA" or target_format == "aa":
            assemble_cmd = str(assemble_cmd) + " -t_aa %s" % targetfile
        elif target_format == "DNA" or target_format == "dna":
            assemble_cmd = str(assemble_cmd) + " -t_dna %s" % targetfile
        for sample in samplenames:
            if sample in file.split('_'):
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
    """In the future, will generate the hybpiper stats command
    (currently uses a static pre-written command for testing)

                Returns
                -------
                stat_cmd : str
                    A string with the run_hybpiper_stats command
    """
    stat_cmd = "hybpiper stats -t_dna test_targets.fasta gene namelist.txt"
    return stat_cmd


def construct_heatmap_command():
    """In the future, will generate the hybpiper generate_heatmap command.
        (currently uses a static pre-written command for testing)

                    Returns
                    -------
                    heatmap_cmd : str
                        A string with the generate_heatmap command
    """
    heatmap_cmd = "hybpiper recovery_heatmap seq_lengths.tsv"
    return heatmap_cmd


def construct_retrieve_commands():
    """In the future, will generate the hybpiper retrieve sequences commands
    and appends them to a list, then returns that list.
    (currently uses static pre-written commands for testing)

                        Returns
                        -------
                        retrieve_cmds : list
                            A list of strings with the hybpiper retrieve
                            sequences commands.
    """
    retrieve_cmds = []
    retrieve_cmd_1 = "hybpiper retrieve_sequences -t_dna test_targets.fasta dna --sample_names namelist.txt --fasta_dir 01_dna_seqs"
    retrieve_cmd_2 = "hybpiper retrieve_sequences -t_dna test_targets.fasta aa --sample_names namelist.txt --fasta_dir 02_aa_seqs"
    retrieve_cmds.append(retrieve_cmd_1)
    retrieve_cmds.append(retrieve_cmd_2)

    return retrieve_cmds


def construct_paralog_command():
    """In the future, will generate the hybpiper paralog retriever command
        (currently uses a static pre-written command for testing)

                    Returns
                    -------
                    paralog_cmd : str
                        A string with the hybpiper paralog retriever command
    """
    paralog_cmd = "hybpiper paralog_retriever namelist.txt -t_dna test_targets.fasta"
    return paralog_cmd


def parseArgvs():
    """This runs the argparse command in order to handle the input arguments
    from the command line.
    """
    parser = argparse.ArgumentParser(prog="generate_namelist.py",
                                     description="Script to create a file "
                                                 "containing the name of "
                                                 "every sample in the readfile ")
    parser.add_argument("-v", "--version", action="version",
                        version="generate_hybpiper_commands.py 1.1.5")
    parser.add_argument("-r", "--readfile", action="store", dest="readfile",
                        help="The location of the input readfile(s)",
                        required=True)
    parser.add_argument("-o", "--output", action="store", dest="output_path",
                        help="The path location where the output file with the new commands should go",
                        required=False)
    parser.add_argument("-t", "--targets", action="store", dest="targetfile",
                        help="The location of the input targetfile",
                        required=True)
    parser.add_argument("-f", "--target_format", action="store",
                        dest="target_format",
                        help="The type of sequences the targets are comprised off, either DNA or Amino acids.",
                        required=True)
    parser.add_argument("-e", "--engine", action="store", dest="search_engine",
                        help="Which method Hybpiper uses to map the sequences to the targets.",
                        required=True)
    parser.add_argument("-i", "--intronerate", action="store",
                        dest="intronerate_bool",
                        help="Whether Hybpiper should extract the introns using intronerate",
                        required=False)
    parser.add_argument("-m", "--heatmap_bool", action="store",
                        dest="heatmap_bool",
                        help="Whether hybpiper should generate a heatmap of the gene recovery",
                        required=False)
    parser.add_argument("-n", "--namelist", action="store",
                        dest="write_namelist",
                        help="Boolean that indicates whether the script should"
                             "write prefixes to a file called 'namelist.txt'",
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
    heatmap_bool = argvs.heatmap_bool
    output_location = argvs.output_path
    write_namelist = argvs.write_namelist

    filenameslist = get_file_names(readfile, False)
    samplenameslist = get_sample_names(filenameslist, True)

    # Generate Hybpiper commands
    cmds = construct_assemble_commands(readfile, filenameslist,
                                       samplenameslist,
                                       target_format, targetfile,
                                       search_method, intronerate_bool)
    #stats_cmd = construct_stats_command()
    #heatmap_cmd = construct_heatmap_command()
    #retrieve_cmds = construct_retrieve_commands()
    #paralog_cmd = construct_paralog_command()

    # write commands to .txt file
    write_commands_to_file(cmds, output_location, overwrite=True)
    # write_commands_to_file([stats_cmd])
    # write_commands_to_file([heatmap_cmd])
    # write_commands_to_file(retrieve_cmds)
    # write_commands_to_file([paralog_cmd])

    # Write namelist.txt
    if write_namelist == 'y' or write_namelist:
        write_commands_to_file(samplenameslist, output_location,
                               "namelist.txt", overwrite=True)


if __name__ == '__main__':
    main()