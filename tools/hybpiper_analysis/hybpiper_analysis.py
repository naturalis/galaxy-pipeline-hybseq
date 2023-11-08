#!/usr/bin/env python

"""
This script needs to be directed to the HybPiper output folder.
It then extracts all the fasta files from the loci that were generated
by the HybPiper pipeline, and sorts them per gene.

These fasta files can then be used to create several MSA's
(multiple sequence alignments) for phylogenetic analysis.
"""
import os
import shutil
import argparse
import subprocess
import sys
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO, AlignIO
from Bio.AlignIO.PhylipIO import PhylipWriter, RelaxedPhylipIterator
from Bio.Application import ApplicationError
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import shlex
import re

global id_length


def concatenate_loci(input_dir, output_dir):
    print("=== concatenate_loci ===")
    print("input_dir:", input_dir)
    print("output_dir:", output_dir)
    """
    Concatenates all the sequences in each locus folder in `input_dir` into a single multi-FASTA file for each locus,
    and saves the files in a new directory called `concatenated` in `output_dir`.

    Parameters
    ----------
    input_dir : str
        The path to the directory containing the locus folders, each of which should contain multiple FNA files with
        sequence data.
    output_dir : str
        The path to the directory where the concatenated files will be saved.

    Returns
    -------
    concat_dir : str
        The path to the `concatenated` directory.

    Raises
    ------
    OSError
        If there is an issue creating the `concatenated` directory.

    Notes
    -----
    This function assumes that the FNA files in each locus folder contain only one sequence each, and that the file
    names are in the format `SAMPLENAME_LOCUS.FNA`, where `SAMPLENAME` is a unique identifier for each sequence and
    `LOCUS` is the name of the locus folder containing the file.
    """
    # Create a folder to store the concatenated files
    concat_dir = os.path.join(output_dir, 'concatenated')
    os.makedirs(concat_dir, exist_ok=True)
    # Iterate over locus folders
    for locus_folder in os.listdir(input_dir):
        locus_folder_path = os.path.join(input_dir, locus_folder)
        if os.path.isdir(locus_folder_path):
            # create a list to store the SeqRecords for each entry in this locus
            seq_records = []
            # iterate over .FNA files in locus folder
            for fna_file in os.listdir(locus_folder_path):
                if fna_file.endswith('.FNA'):
                    fna_file_path = os.path.join(locus_folder_path, fna_file)
                    # get the sample name from the file name
                    sample_name = fna_file.split('_')[0]
                    # iterate over the SeqRecords in the FNA file
                    for seq_record in SeqIO.parse(fna_file_path, 'fasta'):
                        # set the SeqRecord's ID to the sample name
                        seq_record.id = sample_name
                        seq_record.description = ''
                        seq_records.append(seq_record)
            # write the SeqRecords for this locus to the output file
            locus_output_file = os.path.join(concat_dir, locus_folder + '.fa')
            with open(locus_output_file, 'w') as out_file:
                SeqIO.write(seq_records, out_file, 'fasta')
    return concat_dir


def align_loci(input_dir, output_dir):
    global id_length
    print("=== align_loci ===")
    print("input_dir:", input_dir)
    print("output_dir:", output_dir)
    """Performs sequence alignment for each locus in input_dir and 
    saves the resulting alignments in output_dir in the specified format. 
    The sequences will be aligned to the same length, 
    using gaps if necessary.

    Parameters
    ----------
    input_dir : str
        The path to the directory containing the input files.
    output_dir : str
        The path to the directory where the output files will be saved.
    output_format : str, optional
        The format in which to save the output files. Default is '.phy'.

    Returns
    -------
    alignment_dir : str
        The path to the directory where the aligned files were saved.

    """
    # Get the path to the muscle executable
    muscle_path = shutil.which('muscle')
    if muscle_path is None:
        raise ValueError("Muscle executable not found in PATH")

    # Create a subdirectory within the output directory to store the alignment files
    alignment_dir = os.path.join(output_dir, "alignments")
    os.makedirs(alignment_dir, exist_ok=True)

    for concat_file in os.listdir(input_dir):
        input_file = os.path.join(input_dir, concat_file)
        muscle_output_file = os.path.join(
            alignment_dir,
            re.sub("[.](faa|fna|fa)$", ".afa", concat_file)
        )
        phylyp_output_file = os.path.join(
            alignment_dir,
            re.sub("[.](faa|fna|fa)$", ".phy", concat_file)
        )

        # Align sequences using Muscle
        muscle_cline = str(MuscleCommandline(str(muscle_path), input=input_file, out=muscle_output_file))
        # fix muscle command to work with newer version
        muscle_cline = muscle_cline.replace(" -in ", " -align ").replace(" -out ", " -output ")
        print("MUSCLE command:", muscle_cline)
        cmd_tokens = shlex.split(str(muscle_cline))
        proc = subprocess.Popen(cmd_tokens, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        outs, errs = proc.communicate()
        # print("MUSCLE stdout:")
        # print(outs)
        # print("MUSCLE stderr:")
        # print(errs)
        if "invalid" in str(errs):  # FIXME might break if a samplename contains this word
            print("Muscle failed with error")
            exit()
        assert os.path.isfile(cmd_tokens[-1]) and os.stat(cmd_tokens[-1]).st_size != 0, \
            f"MUSCLE output '{cmd_tokens[-1]}' is missing or an empty file"

        # Open the alignment file and align sequences to the same length
        id_length = 0
        with open(muscle_output_file, "r") as alignment_file:
            alignment = AlignIO.read(alignment_file, "fasta")
            alignment_len = alignment.get_alignment_length()
            for record in alignment:
                id_length = max(id_length, len(record.id))
                if len(record.seq) < alignment_len:
                    # If sequence is shorter than the alignment length, add gaps to the end
                    gap_len = alignment_len - len(record.seq)
                    gap_seq = "-" * gap_len
                    record.seq += gap_seq
                elif len(record.seq) > alignment_len:
                    # If sequence is longer than the alignment length, truncate it
                    record.seq = record.seq[:alignment_len]

        # Write the aligned sequences to the same output file in the specified format
        with open(phylyp_output_file, "w") as output_file:
            pw = PhylipWriter(output_file)
            pw.write_alignment(alignment, id_width=id_length)

    return alignment_dir


def merge_alignments(align1, align2):
    """
    Merges two multiple sequence alignments by concatenating sequences
    for matching sequence identifiers. The merged alignment will have gaps
    ('?') for any sequences that are present in only one of the alignments.

    Parameters
    ----------
    align1 : Bio.Align.MultipleSeqAlignment
        The first multiple sequence alignment to merge.
    align2 : Bio.Align.MultipleSeqAlignment
        The second multiple sequence alignment to merge.

    Returns
    -------
    merged_alignment : Bio.Align.MultipleSeqAlignment
        A new multiple sequence alignment containing the merged sequences.
    """
    merged_dict = {}
    for record in align1:
        merged_dict[record.id] = str(record.seq)
    for record in align2:
        if record.id in merged_dict:
            merged_dict[record.id] += str(record.seq)
        else:
            merged_dict[record.id] = "?" * len(align1[0].seq) + str(record.seq)
    for record_id in merged_dict:
        if len(merged_dict[record_id]) < len(align1[0].seq) + len(align2[0].seq):
            merged_dict[record_id] += "?" * (len(align1[0].seq) + len(align2[0].seq) - len(merged_dict[record_id]))
    merged_records = [SeqRecord(Seq(sequence), id=record_id, description="") for record_id, sequence in
                      merged_dict.items()]
    return MultipleSeqAlignment(merged_records)


def create_supermatrix(directory, output_dir, seq_format):
    """
    Reads in multiple phylip-format alignment files in a given directory, concatenates them into a supermatrix alignment
    and writes the concatenated alignment and partition information to output files in the specified output directory.

    Parameters:
    directory (str): the directory containing the input phylip-format alignment files
    output_dir (str): the directory where the output files will be written
    seq_format (str): what format the sequence files are in
    """
    # Initialize the concatenated alignment object and the partition dictionary
    concatenated_alignment = None
    partition_dict = {}
    # Define the starting position for each partition
    start_position = 1
    # Iterate through the files in the directory and process those with the ".phy" extension
    for filename in os.listdir(directory):
        if filename.endswith(".phy"):
            file_path = os.path.join(directory, filename)
            with open(file_path, "r") as fh:
                alignment = list(RelaxedPhylipIterator(fh))[0]  # [0] is to get the full alignment, NOT THE FIRST SEQ
                # it is equivalent to: alignment = AlignIO.read(file_path, "phylip"), but relaxed phylip format
            if concatenated_alignment is None:
                concatenated_alignment = MultipleSeqAlignment(alignment)
            else:
                concatenated_alignment = merge_alignments(concatenated_alignment, alignment)
            # Define the end position for the current partition
            end_position = start_position + len(alignment[0])
            # Add the partition information to the dictionary
            partition_name = "locus{}".format(len(partition_dict) + 1)
            partition_dict[partition_name] = "{}-{}".format(start_position, end_position - 1)
            # Update the starting position for the next partition
            start_position = end_position
    # Define the output file paths
    output_file = os.path.join(output_dir, "concatenated_alignment.phy")
    partition_file = os.path.join(output_dir, "partition.txt")

    # Write the concatenated alignment to the output file
    with open(output_file, "w") as output_handle:
        pw = PhylipWriter(output_handle)
        pw.write_alignment(concatenated_alignment, id_width=id_length)

    # Write the partition information to the partition file
    with open(partition_file, "w") as partition_handle:
        for partition_name, partition_range in partition_dict.items():
            if seq_format == "DNA":
                partition_handle.write("DNA, {} = {}\n".format(partition_name, partition_range))
            else:
                partition_handle.write("AA, {} = {}\n".format(partition_name, partition_range))


def parse_argvs():
    """This runs the argparse command in order to handle the input arguments
    from the command line.
    """
    parser = argparse.ArgumentParser(
        description='Create a phylogenetic tree from enriched sequences of marker genes.')
    parser.add_argument("-v", "--version", action="version",
                        version="hybpiper_analysis.py 1.1.1")
    parser.add_argument('-i', '--input-dir', dest='input_dir', required=True,
                        help='Path to the input directory.')
    parser.add_argument('-o', '--output-dir', dest='output_dir', required=True,
                        help='Path to the output directory.')
    parser.add_argument('-f', '--format', dest='seq_format', required=True,
                        help='What type of sequences? (AA/DNA)')
    argvs = parser.parse_args()
    return argvs


def main():
    argvs = parse_argvs()
    input_dir = argvs.input_dir
    output_dir = argvs.output_dir
    seq_format = argvs.seq_format

    print("Concatenating input files...")
    concatenated_input_dir = concatenate_loci(input_dir, output_dir)
    print("Running MUSCLE...")
    alignment_dir = align_loci(concatenated_input_dir, output_dir)
    print("Generating supermatrix...")
    create_supermatrix(alignment_dir, output_dir, seq_format)
    print("Done!")


if __name__ == '__main__':
    main()
