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
    #muscle_path = os.path.abspath(os.path.join(env_path, "bin", "muscle"))
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
            output_file = concat_file.replace('.FNA', '.afa')
        elif "faa" in concat_file.lower().split("."):
            output_file = concat_file.replace('.FAA', '.afa')
        elif "fa" in concat_file.lower().split("."):
            output_file = concat_file.replace('.fa', '.afa')
        output_path = os.path.join(alignment_dir, output_file)
        muscle_cline = MuscleCommandline(str(muscle_path), input=input_path, out=output_path)
        #print("Running Muscle with command: ", muscle_cline)
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

            # Add a check for file extension before parsing the file as a fasta file
            if not input_path.endswith('.afa'):
                continue

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

    concatenated_input_dir = concatenate_loci(input_dir, output_dir)
    alignment_dir = align_loci(concatenated_input_dir, output_dir)
    supermatrix_file = make_supermatrix(alignment_dir, output_dir)
    #make_partition_file(supermatrix_file, output_dir)

    #print("Input dir: " + str(input_dir))
    #print("Output dir: " + str(output_dir))
    #print("Alignment_dir: " + str(alignment_dir))
    # print("Supermatrix_file: " + str(supermatrix_file))


if __name__ == '__main__':
    main()