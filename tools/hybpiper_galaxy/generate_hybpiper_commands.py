#!/usr/bin/env python

"""
This script takes several inputs using the argparse command, with which it
generates the commands needed to run Hybpiper.

It works in conjunction with- and is called by the hybpiper_galaxy.sh shell
script to create a file with the complete commands.

The shell script would then go over each command in that text file one by one
and execute them to run the Hybpiper pipeline as intended.
"""
import os
import argparse
import re
import zipfile
import shutil


def get_filenames(zip_files_list: list[str]) -> tuple[list[str], list[str]]:
    """
    :param zip_files_list: list of files in the zip file
    :return: the list of filenames (without paths) + filenames with path in the zip file
    """
    full_filenames_list = [n for n in zip_files_list
                           if (n.endswith(".fastq") or n.endswith(".fastq.gz"))]
    assert full_filenames_list, "Fastq file cannot be found in zip file"
    # check that all files are located same directory or at the root
    store_dir = set([os.path.dirname(n) for n in full_filenames_list])
    assert len(store_dir) == 1, "Fastq files found in different location inside the zip file"
    filenames_list = [os.path.basename(n) for n in full_filenames_list]

    R1_n = [n for n in filenames_list if "_R1_" in n]
    R2_n = [n for n in filenames_list if "_R2_" in n]
    assert R1_n and R2_n, 'All read files must contain "_R1_" or "_R2_" in their names'
    # assert len(R1_n) == len(R2_n), 'Different number of files containing "_R1_" and "_R2_" in read dir'
    for n in R1_n:
        assert n.replace("_R1_", "_R2_") in R2_n, f'Cannot find corresponding R2 file of "{n}"'
    for n in R2_n:
        assert n.replace("_R2_", "_R1_") in R1_n, f'Cannot find corresponding R1 file of "{n}"'

    return filenames_list, full_filenames_list


def extract_read_files(read_zip_file: str,
                       output_dir: str) -> list[str]:
    assert zipfile.is_zipfile(read_zip_file), "Input Zip file is not a valid ZIP file"
    with zipfile.ZipFile(read_zip_file, 'r') as ziph:
        names = ziph.namelist()
        print("Names in zip file:", len(names))
        filenames, full_names = get_filenames(names)
        print("FASTQ Filenames:", len(filenames))
        # only extract the fastq found
        ziph.extractall(path=output_dir, members=full_names, pwd=None)
        # move the fastq files at the root of the output_dir
        for root, dirs, files in os.walk(output_dir):
            for f in files:
                if os.path.join(root, f) != os.path.join(output_dir, f):
                    shutil.move(os.path.join(root, f), output_dir)

    return filenames


def write_commands_to_file(commands: list[str],
                           output_dir: str,
                           output_name: str,
                           append: bool = False) -> str:
    """Iterates through a list, and writes each element to a new line in the
    specified output file.

                Parameters
                ----------
                commands : list
                    A list with all the elements the function needs to write
                    to the output file.
                output_dir : str
                    A string with the path to where the output file should be
                    by default the output directory is the current
                    working directory
                output_name : str
                    by default the output file is a txt file called
                    'cmdfile.txt' in the same directory as the script.
                append : bool
                    A boolean indicating whether the outputfile should be
                    reused. If True, output is appended to the output file, if
                    False, the output file is emptied before being written to.
    """
    output_file = os.path.join(output_dir, output_name)
    if not append:
        with open(output_file, "w") as fw:
            fw.close()
    with open(output_file, "a") as outfile:
        outfile.write("\n".join(commands) + "\n")
    return output_file


def construct_assemble_commands(read_dir: str,
                                output_dir: str,
                                filenames: list[str],
                                samplenames: list[str],
                                target_format: str,
                                targetfile: str,
                                search_engine: str,
                                no_intronerate: bool,
                                timeout: int) -> list[str]:
    """The function constructs the hybpiper assemble command
    It does this by appending a template string with the proper flags according
    to the state of the arguments. Using this method, it assembles one
    command for each sample for every read file. Then returns a list
    of all the generated commands.
        Parameters
        ----------
        read_dir : str
            A string resembling the path to the folder with the read
            fastq files.
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
            A string corresponding to the type of method hybpiper
            should use for mapping reads.
        no_intronerate : str
            Boolean to add the flag --no_intronerate
            to the hybpiper assemble command.
        timeout : int
            Represents an integer x, that allows hybpiper to kill processes that take x percent longer than the average,
            preventing the tool from getting stuck.

        Returns
        -------
        assemble_cmds : list
            a list of the generated hybpiper assemble commands in the
            form of strings.
    """
    assemble_cmds = []
    for file in filenames:
        assemble_cmd = "hybpiper assemble -r %s -o %s" % (os.path.join(read_dir, file), output_dir)

        if target_format == "AA":
            assemble_cmd += " -t_aa %s" % targetfile
        else:  # dna
            assemble_cmd += " -t_dna %s" % targetfile

        for sample in samplenames:
            assemble_cmd += " --prefix %s" % sample

        if search_engine != "blastx":
            assemble_cmd += " --%s" % search_engine

        if no_intronerate:
            assemble_cmd += " --no_intronerate"

        if timeout:
            assemble_cmd += " --timeout_assemble %s" % timeout
        # assemble_cmd = str(assemble_cmd) + " --hybpiper_dir %s" % str(output_location)
        assemble_cmds.append(assemble_cmd)
    return assemble_cmds


def construct_stats_command() -> str:
    """In the future, will generate the hybpiper stats command
    (currently uses a static pre-written command for testing)

                Returns
                -------
                stat_cmd : str
                    A string with the run_hybpiper_stats command
    """
    return "hybpiper stats -t_dna test_targets.fasta gene namelist.txt"


def construct_heatmap_command() -> str:
    """In the future, will generate the hybpiper generate_heatmap command.
        (currently uses a static pre-written command for testing)

                    Returns
                    -------
                    heatmap_cmd : str
                        A string with the generate_heatmap command
    """
    return "hybpiper recovery_heatmap seq_lengths.tsv"


def construct_retrieve_commands() -> list[str]:
    """In the future, will generate the hybpiper retrieve sequences commands
    and appends them to a list, then returns that list.
    (currently uses static pre-written commands for testing)

                        Returns
                        -------
                        retrieve_cmds : list
                            A list of strings with the hybpiper retrieve
                            sequences commands.
    """
    # FIXME Hard coded filenames
    retrieve_cmds = [
        "hybpiper retrieve_sequences -t_dna test_targets.fasta dna --sample_names namelist.txt --fasta_dir 01_dna_seqs",
        "hybpiper retrieve_sequences -t_dna test_targets.fasta aa --sample_names namelist.txt --fasta_dir 02_aa_seqs"
    ]

    return retrieve_cmds


def construct_paralog_command() -> str:
    """In the future, will generate the hybpiper paralog retriever command
        (currently uses a static pre-written command for testing)

                    Returns
                    -------
                    paralog_cmd : str
                        A string with the hybpiper paralog retriever command
    """
    return "hybpiper paralog_retriever namelist.txt -t_dna test_targets.fasta"


def parse_argvs() -> argparse.Namespace:
    """This runs the argparse command in order to handle the input arguments
    from the command line.
    """
    parser = argparse.ArgumentParser(prog="generate_namelist.py",
                                     description="Script to create a file "
                                                 "containing the name of "
                                                 "every sample in the readfile ")
    parser.add_argument("-v", "--version",
                        action="version",
                        version="generate_hybpiper_commands.py 1.1.9")
    parser.add_argument("--read_zip",
                        help="The location of the input zip file containing the reads (fastq files)",
                        required=True)
    parser.add_argument("--read_dir",
                        help="The location where to extract the read file(s)",
                        required=True)
    parser.add_argument("--output_dir",
                        help="The path location where the output file with the new commands should go",
                        required=False)
    parser.add_argument("--targets",
                        dest="targetfile",
                        help="The location of the input targetfile",
                        required=True)
    parser.add_argument("--target_format",
                        choices=["DNA", "AA"],
                        help="The type of sequences the targets are comprised off, either DNA or Amino acids.",
                        required=True)
    parser.add_argument("--engine",
                        dest="search_engine",
                        choices=["blastx", "bwa", "diamond"],
                        help="Which method Hybpiper uses to map the sequences to the targets.",
                        required=True)
    parser.add_argument("--no_intronerate",
                        action="store_true",
                        help="Whether Hybpiper should not extract the introns using intronerate"
                             " (it is performed by Hybpiper by default)",
                        required=False)
    parser.add_argument("--heatmap",
                        action="store_true",
                        help="Whether hybpiper should generate a heatmap of the gene recovery",
                        required=False)
    parser.add_argument("--timeout",
                        type=int,
                        help="Enter a whole number X, the program will kill processes that "
                             "take X percent longer than average. Use this if jobs get stuck. Leave empty for default.",
                        required=False)
    parser.add_argument("--write_namelist",
                        action="store_true",
                        help="Boolean that indicates whether the script should"
                             "write prefixes to a file called 'namelist.txt'",
                        required=False)
    argvs = parser.parse_args()
    return argvs


def main():
    argvs = parse_argvs()
    zip_read_file = argvs.read_zip
    read_dir = argvs.read_dir
    targetfile = argvs.targetfile
    target_format = argvs.target_format
    search_method = argvs.search_engine
    no_intronerate = argvs.no_intronerate
    heatmap = argvs.heatmap
    output_dir = argvs.output_dir
    write_namelist = argvs.write_namelist
    timeout = argvs.timeout

    filenames = extract_read_files(zip_read_file, read_dir)
    samplenames = list(set([re.split("_R[12]_", n)[0] for n in filenames]))
    assert len(samplenames) == len(filenames) / 2  # in principle, it must be the case

    # Generate Hybpiper commands
    cmds = construct_assemble_commands(read_dir,
                                       os.path.join(output_dir, "assemble"),
                                       filenames,
                                       samplenames,
                                       target_format,
                                       targetfile,
                                       search_method,
                                       no_intronerate,
                                       timeout)

    # stats_cmd = construct_stats_command()
    # heatmap_cmd = construct_heatmap_command()
    # retrieve_cmds = construct_retrieve_commands()
    # paralog_cmd = construct_paralog_command()

    # write commands to .txt file
    cmd_file = write_commands_to_file(cmds,
                                      output_dir,
                                      "cmdfile.txt",
                                      append=False)
    assert os.stat(cmd_file).st_size != 0, "Commands file is empty"

    # write_commands_to_file([stats_cmd])
    # write_commands_to_file([heatmap_cmd])
    # write_commands_to_file(retrieve_cmds)
    # write_commands_to_file([paralog_cmd])

    # Write namelist.txt
    if write_namelist:
        samplename_file = write_commands_to_file(samplenames,
                                                 output_dir,
                                                 "namelist.txt",
                                                 append=False)
        assert os.stat(samplename_file).st_size != 0, "Sample names file is empty"


if __name__ == '__main__':
    main()
