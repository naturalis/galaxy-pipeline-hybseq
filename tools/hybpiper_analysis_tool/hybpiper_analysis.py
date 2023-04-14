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
from Bio.Application import ApplicationError


def concatenate_loci(input_dir, output_dir):
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


def align_loci(input_dir, output_dir, output_format='.phy'):
    """Performs sequence alignment for each locus in input_dir and
    saves the resulting alignments in output_dir in the specified format.
    The sequences will be aligned to the same length, using gaps if necessary.
    """
    # Get the path to the muscle executable
    muscle_path = shutil.which('muscle')
    if muscle_path is None:
        raise ValueError("Muscle executable not found in PATH")

    # Create a subdirectory within the output directory to store the alignment files
    alignment_dir = os.path.join(output_dir, "alignments")
    os.makedirs(alignment_dir, exist_ok=True)

    concat_dir = os.path.join(output_dir, "concatenated")
    concat_path = os.listdir(concat_dir)

    for concat_file in concat_path:
        input_path = os.path.join(concat_dir, concat_file)
        if "fna" in concat_file.lower().split("."):
            output_file = concat_file.replace('.fna', output_format)
        elif "faa" in concat_file.lower().split("."):
            output_file = concat_file.replace('.faa', output_format)
        elif "fa" in concat_file.lower().split("."):
            output_file = concat_file.replace('.fa', output_format)
        else:
            continue

        output_path = os.path.join(alignment_dir, output_file)
        # Align sequences using Muscle
        muscle_cline = MuscleCommandline(str(muscle_path), input=input_path, out=output_path)
        try:
            stdout, stderr = muscle_cline()
        except ApplicationError as e:
            print("Muscle failed with error:", e.stderr)
            continue

        # Open the alignment file and align sequences to the same length
        with open(output_path, "r") as alignment_file:
            alignment = AlignIO.read(alignment_file, "fasta")
            alignment_len = alignment.get_alignment_length()
            for record in alignment:
                if len(record.seq) < alignment_len:
                    # If sequence is shorter than the alignment length, add gaps to the end
                    gap_len = alignment_len - len(record.seq)
                    gap_seq = "-" * gap_len
                    record.seq += gap_seq
                elif len(record.seq) > alignment_len:
                    # If sequence is longer than the alignment length, truncate it
                    record.seq = record.seq[:alignment_len]

        # Write the aligned sequences to the same output file in the specified format
        with open(output_path, "w") as output_file:
            AlignIO.write(alignment, output_file, "phylip")

    return alignment_dir


import os

def create_supermatrix(input_folder, output_folder):
    # Create a list of all the input file paths
    input_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.phy')]

    # Determine the number of loci and the length of the alignment for each file
    num_loci = []
    alignment_lengths = []
    for file in input_files:
        with open(file, 'r') as f:
            num_sequences, alignment_length = f.readline().split()
            num_sequences = int(num_sequences)
            alignment_length = int(alignment_length)
            num_missing = 0
            for line in f:
                if line.startswith(' '):
                    num_missing += 1
            num_samples = num_sequences - num_missing
            num_loci.append(num_samples)
            alignment_lengths.append(alignment_length)

    # Create the partition definition block
    num_partitions = len(num_loci)
    partition_sizes = ['{}-{}'.format(sum(num_loci[:i])+1, sum(num_loci[:i+1])) for i in range(num_partitions)]
    partition_types = ['DNA' for i in range(num_partitions)]
    partition_def = '{} {} {}\n'.format(num_partitions, ' '.join(partition_types), ' '.join(partition_sizes))

    # Create the partition.txt file
    with open(os.path.join(output_folder, 'partition.txt'), 'w') as f:
        for i in range(num_partitions):
            f.write('{} locus{}\n'.format(partition_sizes[i], i+1))

    # Create the output supermatrix file
    with open(os.path.join(output_folder, 'supermatrix.phy'), 'w') as out_file:
        out_file.write(partition_def)
        for i in range(num_partitions):
            with open(input_files[i], 'r') as in_file:
                in_file.readline()
                for line in in_file:
                    if line.startswith(' '):
                        out_file.write('?' * alignment_lengths[i] + '\n')
                    else:
                        out_file.write(line)



def parse_argvs():
    """This runs the argparse command in order to handle the input arguments
    from the command line.
    """
    parser = argparse.ArgumentParser(
        description='Create a phylogenetic tree from enriched sequences of marker genes.')
    parser.add_argument('-i', '--input-dir', dest='input_dir', required=True,
                        help='Path to the input directory.')
    parser.add_argument('-o', '--output-dir', dest='output_dir', required=True,
                        help='Path to the output directory.')
    argvs = parser.parse_args()
    return argvs


def main():
    LINE_CLEAR = '\x1b[2K'
    argvs = parse_argvs()
    input_dir = argvs.input_dir
    output_dir = argvs.output_dir

    print("Concatenating input files...", end="\r")
    concatenated_input_dir = concatenate_loci(input_dir, output_dir)
    print(end=LINE_CLEAR)
    print("Running MUSCLE...", end="\r")
    alignment_dir = align_loci(concatenated_input_dir, output_dir)
    print(end=LINE_CLEAR)
    print("Generating supermatrix...", end="\r")
    create_supermatrix(alignment_dir, output_dir)
    print(end=LINE_CLEAR)
    print("Done!", end="\r")


if __name__ == '__main__':
    main()