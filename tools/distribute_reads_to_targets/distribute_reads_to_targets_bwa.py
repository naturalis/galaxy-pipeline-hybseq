#!/usr/bin/env python

'''
Edited by Jeremy van Veen commissioned by Naturalis for use in Galaxy
Original Author: https://github.com/mossmatters
Version 1.4.8

Usage:
-------------------
python3 distribute_reads_to_targets_bwa.py -b path_to_bamfile -r
raw_reads_folder_path

python3 distribute_reads_to_targets_bwa.py -b Galaxy82-[Merged_BAM].bam -r backup

Description:
----------------
This script is part of a pipeline to extract phylogenetically-useful
sequences from Illumina data using the targeted (liquid-phase) sequence
enrichment approach.

After a BWA search of the raw reads against the target sequences,
the reads need to be sorted according to the successful hits.
This script takes the BWA output (BAM format) and the raw read files,
and distributes the reads into FASTA files ready for assembly.

If there are multiple results (for example, one for each read direction),
concatenate them prior to sorting.

Notes:
    Currently only tested with paired data, implementation for unpaired
    requires a slight amount of modification in the future.
'''

import sys, os, errno, subprocess, re, argparse
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def getFilenames(directory):
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
            file = file.replace('/', "\\")
            filenames_list = file.split("\\")
            filenames.append(filenames_list[-1:][0])
    return filenames


def mkdir_p(path):
    """Attempts to create output directories if possible.

            Parameters
            ----------
            path : str
                The path to the output directory
            """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def read_sorting(bamfilename):
    """uses samtools to read the BAMfile and sort the reads to the targets
    then updates a dictionary for each hit.

        Parameters
        ----------
        bamfilename : str
            a string with the name of the bamfile

        Returns
        -------
        read_hit_dict : dictionary
            a dictionary containing the read hits towards the target
            sequences in the BAM file
        """
    samtools_cmd = "samtools view -F 4 {}".format(bamfilename)
    child = subprocess.Popen(samtools_cmd, shell=True, stdout=subprocess.PIPE,
                             universal_newlines=True)
    bwa_results = child.stdout.readlines()

    read_hit_dict = {}
    for line in bwa_results:
        line = line.split()
        readID = line[0]
        target = line[2].split("-")[-1]
        if readID in read_hit_dict:
            if target not in read_hit_dict[readID]:
                read_hit_dict[readID].append(target)
        else:
            read_hit_dict[readID] = [target]
    return read_hit_dict


def write_paired_seqs(target, ID1, Seq1, ID2, Seq2, single=True):
    """Writes results towards the output fasta files

            Parameters
            ----------
            target : str
                a string containing the output path
            ID1 : str
                a string containing ID 1
            ID2 : str
                a string containing ID 2
            Seq1 : str
                a string containing sequence1
            Seq2 : str
                a string containing sequence2
            single : bool, optional
                boolean to specify whether the reads are single-end or paired
                (default=True)
            """
    path = "fastafiles/"
    mkdir_p(path)
    if single:
        outfile = open(
            os.path.join(path, "{}_interleaved.fasta".format(target)), 'a')
        outfile.write(">{}\n{}\n".format(ID1, Seq1))
        outfile.write(">{}\n{}\n".format(ID2, Seq2))
        outfile.close()
    else:
        outfile1 = open(os.path.join(path, "{}_1.fasta".format(target)), 'a')
        outfile1.write(">{}\n{}\n".format(ID1, Seq1))
        outfile2 = open(os.path.join(path, "{}_2.fasta".format(target)), 'a')
        outfile2.write(">{}\n{}\n".format(ID2, Seq2))
        outfile1.close()
        outfile2.close()


def write_single_seqs(target, ID1, Seq1):
    """Distributing targets from single-end sequencing

               Parameters
               ----------
               target : str
                   a string containing the output path
               ID1 : str
                   a string containing ID 1
               Seq1 : str
                   a string containing sequence1
               """
    path = "fastafiles/"
    mkdir_p(path)
    outfile = open(os.path.join(path, "{}_unpaired.fasta".format(target)),
                   'a')
    outfile.write(">{}\n{}\n".format(ID1, Seq1))
    outfile.close()


def distribute_reads(readfile, pair_list, read_hit_dict, single=True):
    """uses samtools to read the BAMfile and sort the reads to the targets
        then updates a dictionary for each hit.

            Parameters
            ----------
            readfile : str
                path to the read folder
            read_hit_dict : dictionary
                a dictionary containing the hits from the target function
            singe : bool, optional
                a boolean indicating whether the reads are single-end
                or paired (default=True)
            """
    filenames = pair_list
    num_reads_to_write = len(read_hit_dict)
    if num_reads_to_write != 0:
        iterator1 = FastqGeneralIterator(open("".join([readfile, '/',
                                                       filenames[0]])))
        reads_written = 0
        sys.stderr.write("Read distributing progress:\n")

        if len(filenames) == 1:
            for ID1_long, Seq1, Qual1 in iterator1:
                ID1 = ID1_long.split()[0]
                if ID1.endswith("\1") or ID1.endswith("\2"):
                    ID1 = ID1[:-2]
                if ID1 in read_hit_dict:
                    for target in read_hit_dict[ID1]:
                        write_single_seqs(target, ID1, Seq1)
                        reads_written += 1
                j = (reads_written + 1) / num_reads_to_write
                if int(100 * j) % 5 == 0:
                    sys.stderr.write("\r")
                    sys.stderr.write("[%-20s] %d%%" % ('=' * int(20 * j),
                                                       100 * j))
                    sys.stderr.flush()
            sys.stderr.write("\n")
            return

        elif len(filenames) == 2:
            iterator2 = FastqGeneralIterator(open("".join([readfile, '/',
                                                           filenames[1]])))

        for ID1_long, Seq1, Qual1 in iterator1:

            ID2_long, Seq2, Qual2 = next(iterator2)

            ID1 = ID1_long.split()[0]
            if ID1.endswith("/1") or ID1.endswith("/2"):
                ID1 = ID1[:-2]

            ID2 = ID2_long.split()[0]
            if ID2.endswith("/1") or ID2.endswith("/2"):
                ID2 = ID2[:-2]

            if ID1 in read_hit_dict:
                for target in read_hit_dict[ID1]:
                    write_paired_seqs(target, ID1, Seq1, ID2, Seq2)
                reads_written += 1
            elif ID2 in read_hit_dict:
                for target in read_hit_dict[ID2]:
                    write_paired_seqs(target, ID1, Seq1, ID2, Seq2)
                reads_written += 1
            j = (reads_written + 1) / num_reads_to_write
            if int(100 * j) % 5 == 0:
                sys.stderr.write("\r")
                sys.stderr.write("[%-20s] %d%%" % ('=' * int(20 * j), 100 * j))
                sys.stderr.flush()
        sys.stderr.write("\n")

# The argvs function.
def parseArgvs():
    parser = argparse.ArgumentParser(prog="distribute_reads_"
                                          "to_targets_bwa.py",
                                     description="Script to distribute "
                                                 "reads from a BAM file "
                                                 "and the raw read files "
                                                 "and  distributes them "
                                                 "into fasta files for"
                                                 "assembly")
    parser.add_argument("-v", "--version", action="version",
                        version="distribute_reads_to_targets_bwa 1.4.8")
    parser.add_argument("-b", "--bamfile", action="store", dest="bamfile",
                        help="The location of the input bamfile",
                        required=True)
    parser.add_argument("-r", "--readfile", action="store", dest="readfile",
                        help="The location of the input readfile(s)",
                        required=True)
    argvs = parser.parse_args()
    return argvs

## Scrapped Function
# def getSampleNames(filenames):
#     samplenames_list, pair_dict = [], {}
#     for filename in filenames:
#         part_index = 0
#         filename_parts = filename.split("_")
#         if "R1" in filename_parts:
#             part_index = filename_parts.index("R1")
#         elif "R2" in filename_parts:
#             part_index = filename_parts.index("R2")
#         samplename = filename_parts[:part_index][0]
#         samplenames_list.append(samplename)


def makePairs(filenames):
    """Receives the list of filenames created by getFilenames and pairs them
    together in a list made of tuples.

           Parameters
           ----------
           filenames : list
               list with all filenames from the input read directory as
               strings

           Returns
           -------
           pair_list : list
               a list containing pairs of readfile names as tuples.
           None
               in case of uneven number of reads in input file
           """
    if len(filenames) % 2 == 0:
        pair_list = []
        filenames_forward = filenames[0::2]
        filenames_reverse = filenames[1::2]
        for counter in range(len(filenames_forward)):
            pair = (filenames_forward[counter], filenames_reverse[counter])
            pair_list.append(pair)
        return pair_list
    else:
        print("File contains an uneven number of reads, "
              "unable to create pairs. Return None")
        return None


def main():
    argvs = parseArgvs()
    bamfilename = argvs.bamfile
    read_zip_file = argvs.readfile
    filenames = getFilenames(read_zip_file)
    pair_list = makePairs(filenames)

    read_hit_dict = read_sorting(bamfilename)
    print("Unique reads with hits: {}".format(len(read_hit_dict)))
    for pair in pair_list:
        distribute_reads(read_zip_file, pair, read_hit_dict, single=True)


if __name__ == "__main__": main()

