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
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Application import ApplicationError


def align_loci(input_dir, output_dir):

    """Performs sequence alignment for each locus in input_dir and
    saves the resulting alignments in output_dir.

    """

    ### Get path to muscle.exe ###

    # Get the name of the current Conda environment
    env_name = os.environ.get("CONDA_DEFAULT_ENV")
    # Run the conda env list command to get a list of Conda environments and their paths
    output = subprocess.check_output([sys.executable, "-m", "conda", "env", "list"])
    # Convert the byte string output to a regular string and split it into lines
    output_lines = output.decode().split("\n")
    # Find the line that contains the current environment name
    env_line = next(line for line in output_lines if env_name in line)
    # Extract the path to the environment from the environment line
    env_path = env_line.split()[1]
    # Construct the path to the muscle executable
    muscle_path = shutil.which('muscle')
    if muscle_path is None:
        raise ValueError("Muscle executable not found in PATH")
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    # Create a subdirectory within the output directory to store the alignment files
    alignment_dir = os.path.join(output_dir, "alignments")
    os.makedirs(alignment_dir, exist_ok=True)
    # For each locus folder, execute the other code
    for locus_dir in os.listdir(input_dir):
        locus_path = os.path.join(input_dir, locus_dir)
        if not os.path.isdir(locus_path):
            continue
        # Check if file extention is among allowed file types.
        fasta_files = [f for f in os.listdir(locus_path) if f.endswith('.FNA') or f.endswith('.FAA')]
        if not fasta_files:
            continue
        # Create a subdirectory within the alignment directory for this locus
        locus_alignment_dir = os.path.join(alignment_dir, locus_dir)
        os.makedirs(locus_alignment_dir, exist_ok=True)
        converted_files = []
        for fasta_file in fasta_files:
            input_path = os.path.join(locus_path, fasta_file)
            if "fna" in fasta_file.lower().split("."):
                output_file = fasta_file.replace('.FNA', '.afa')
            elif "fa" in fasta_file.lower().split("."):
                output_file = fasta_file.replace('.fa', '.afa')
            output_path = os.path.join(locus_alignment_dir, output_file)
            muscle_cline = MuscleCommandline(str(muscle_path), input=input_path, out=output_path)
            print("Running Muscle with command: ", muscle_cline)
            try:
                stdout, stderr = muscle_cline()
            except ApplicationError as e:
                print("Muscle failed with error:", e.stderr)
                continue
    return alignment_dir


def make_supermatrix(alignment_dir, output_dir):
    """Concatenates alignment files from alignment_dir into a
    supermatrix and saves the resulting alignment in output_dir.
    """
    alignments = []
    for locus_dir in os.listdir(alignment_dir):
        locus_path = os.path.join(alignment_dir, locus_dir)
        if not os.path.isdir(locus_path):
            continue
        alignment_files = [f for f in os.listdir(locus_path) if f.endswith('.afa')]
        if not alignment_files:
            continue
        for alignment_file in alignment_files:
            input_path = os.path.join(locus_path, alignment_file)
            with open(input_path) as f:
                alignment = AlignIO.read(f, 'fasta')
                alignments.append(alignment)
    supermatrix = MultipleSeqAlignment(alignments)
    for record in supermatrix:
        if 'X' in record.seq:
            record.seq = SeqRecord(seq=record.seq.ungap('X'), id=record.id, name=record.name, description=record.description).seq
    supermatrix_file = os.path.join(output_dir, 'supermatrix.phy')
    with open(supermatrix_file, 'w') as f:
        SeqIO.write(supermatrix, f, 'phylip')
    return supermatrix_file


def make_partition_file(supermatrix_file, output_dir):
    """Generates a partition file for supermatrix_file and saves it
    in output_dir.
    """
    with open(supermatrix_file) as f:
        alignment = AlignIO.read(f, 'phylip')
    partition_file = os.path.join(output_dir, 'partitions.txt')
    with open(partition_file, 'w') as f:
        for i, column in enumerate(alignment):
            start = i+1
            end = i+2
            gene_name = f'gene{i+1}'
            f.write(f'{start}-{end} {gene_name}\n')


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
    argvs = parse_argvs()
    input_dir = argvs.input_dir
    output_dir = argvs.output_dir

    alignment_dir = align_loci(input_dir, output_dir)
    #supermatrix_file = make_supermatrix(alignment_dir, output_dir)
    #make_partition_file(supermatrix_file, output_dir)

    #print("Input dir: " + str(input_dir))
    #print("Output dir: " + str(output_dir))
    #print("Alignment_dir: " + str(alignment_dir))
    # print("Supermatrix_file: " + str(supermatrix_file))


if __name__ == '__main__':
    main()