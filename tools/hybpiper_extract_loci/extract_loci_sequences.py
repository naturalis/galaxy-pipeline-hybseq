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
import re


def extract_loci(input_dir, output_dir, file_option="FNA"):
    """Goes through the HybPiper output files and searches for the specified
    file extentions, then copies every sequence file that matches the search
    extention to a seperate folder before calling the sort_files() function.
        Parameters
        ----------
        input_dir : str
            A string with the specified absolute path to the folder
            where the hybpiper output directory can be found
        output_dir : str
            A string with the specified absolute path to the directory
            where the sequence files should be extracted to.
        file_option : str
            A string indicating what extention to look for options are
            ('FNA', 'FAA' and 'intron')
    """
    allowed_options = ["FNA", "FAA", "intron"]
    fna_dir_ext, path_dir_ext, output_file_dir = file_option, file_option, \
        file_option

    if file_option in allowed_options:
        if file_option == "intron":
            fna_dir_ext = file_option
            path_dir_ext = "fasta"
            output_file_dir = "fasta"
        else:
            pass
    else:
        print("The -t flag had an illegal option, allowed options: %s"
              % allowed_options)
        exit()
    # Create the output directory if it doesn't already exist
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # Loop through each sample directory
    for sample_dir in os.listdir(input_dir):
        # Check if the directory is a sample directory
        if not os.path.isdir(os.path.join(input_dir, sample_dir)):
            continue
        # Loop through each gene_xxx folder within the sample directory
        for gene_dir in os.listdir(os.path.join(input_dir, sample_dir)):
            if not os.path.isdir(os.path.join(input_dir, sample_dir, gene_dir)):
                continue
            # Navigate to the FNA/FAA/intron folder within the gene_xxx folder
            fna_dir = os.path.join(input_dir, sample_dir, gene_dir, sample_dir,
                                   "sequences", "%s" % fna_dir_ext)
            if os.path.isdir(fna_dir):
                # Loop through each file within the FNA folder
                for fna_file in os.listdir(fna_dir):
                    # Check if the file is of the specified type
                    if not fna_file.endswith(".%s" % path_dir_ext):
                        continue
                    # Create the output path for the file
                    if file_option == "intron":
                        output_file = os.path.join(output_dir, fna_file)
                    else:
                        output_file = os.path.join(output_dir,
                                                   sample_dir + "_" + gene_dir
                                                   + "." + "%s"
                                                   % output_file_dir)
                    #print(output_file)
                    # Copy the file to the output directory
                    shutil.copy(os.path.join(fna_dir, fna_file), output_file)
                    # Sort files by gene
                    sort_files(output_dir)


def sort_files(loci_folder_path):
    """Iterates through the folder specified in the arguments and sorts the
    sequence files that are present by the gene number in the filename.
    A folder is created for each gene number present and numbers are then
    sorted into said folders.
            Parameters
            ----------
            loci_folder_path : str
                A string with the specified absolute path to the folder
                where the sequence files are located, that need to be
                sorted per gene.
    """
    # Create a generator expression to loop over
    # all the loci files one at a time
    letter_pattern = r'[a-zA-Z]'
    loci_files = (f for f in os.listdir(loci_folder_path) if
                  os.path.isfile(os.path.join(loci_folder_path, f)))
    # Loop over all the loci files
    for loci_file in loci_files:
        # Remove extention
        file_parts = loci_file.split(".")[0]
        # Split the file name into its components
        file_parts = file_parts.split('_')
        sample_name = file_parts[0]
        # Remove any and all letters from the gene_number extention
        gene_number = re.sub(letter_pattern, '', file_parts[1])
        # Create the folder name
        gene_folder_name = 'gene{}_{}'.format(gene_number,
                                              loci_file.split(".")[1])
        # Create the gene folder if it doesn't already exist
        gene_folder_path = os.path.join(loci_folder_path, gene_folder_name)
        os.makedirs(gene_folder_path, exist_ok=True)
        # Move the file into the gene folder
        old_file_path = os.path.join(loci_folder_path, loci_file)
        shutil.move(old_file_path, gene_folder_path)


def count_genes_per_sample(loci_folder_path, output_file_path):
    # Create a dictionary to hold the counts for each sample
    sample_counts = {}

    # Loop over all the gene folders
    for gene_folder_name in os.listdir(loci_folder_path):
        # Ignore any files that are not directories
        gene_folder_path = os.path.join(loci_folder_path, gene_folder_name)
        if not os.path.isdir(gene_folder_path):
            continue

        # Loop over all the files in the gene folder
        for file_name in os.listdir(gene_folder_path):
            # Split the file name into its components
            file_parts = file_name.split('_')
            sample_name = file_parts[0]

            # Increment the count for this sample
            if sample_name in sample_counts:
                sample_counts[sample_name] += 1
            else:
                sample_counts[sample_name] = 1

    # Write the results to a text file
    with open(str(output_file_path) + "/statistics.txt", 'w') as output_file:
        output_file.write('Amount of sequence files per sample\n')
        output_file.write('---\n')
        for sample_name, count in sample_counts.items():
            output_file.write('{}:\t{}\n'.format(sample_name, count))


def parse_argvs():
    """This runs the argparse command in order to handle the input arguments
    from the command line.
    """
    parser = argparse.ArgumentParser(prog="extract_loci_sequences.py",
                                     description="Script to search the loci "
                                                 "sequences from standard "
                                                 " HybPiper output folder.")
    parser.add_argument("-v", "--version", action="version",
                        version="extract_loci_sequences.py 1.0.7")
    parser.add_argument("-f", "--hybpiper_folder", action="store",
                        dest="input_path",
                        help="The path/location where the script should "
                             "extract the loci from",
                        required=True)
    parser.add_argument("-o", "--output_folder", action="store",
                        dest="output_path",
                        help="The path/location where the script should "
                             "extract the loci to",
                        required=True)
    parser.add_argument("-t", "--filetype", action="store", dest="file_option",
                        help="options: FAA, FNA, intron the type of sequence "
                             "file the script is supposed to extract",
                        default="FNA")
    parser.add_argument("-c", "--count", action="store", dest="count_bool",
                        help="Whether or not to count the amount of sequence"
                             "files per sample (may take longer if True)",
                        default="False")
    argvs = parser.parse_args()
    return argvs


def main():
    argvs = parse_argvs()
    input_path = argvs.input_path
    output_path = argvs.output_path
    file_option = argvs.file_option
    count = argvs.count_bool

    extract_loci(input_path, output_path, file_option)

    if count == "True":
        count_genes_per_sample(output_path, output_path)


if __name__ == '__main__':
    main()
