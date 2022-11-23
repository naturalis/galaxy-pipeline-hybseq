#!/usr/bin/env python

'''
Edited by Jeremy van Veen for Naturalis for use in Galaxy
Version 1.3.2

Use: python3 distribute_reads_to_targets_bwa.py
-b Galaxy82-[Merged_BAM].bam -r test_dataset/test_reads.fastq/EG30_test
'''

import sys, os, errno, subprocess, re, argparse
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

"""
This script is part of a pipeline to extract phylogenetically-useful sequences from 
Illumina data using the targeted (liquid-phase) sequence enrichment approach.

After a BWA search of the raw reads against the target sequences, the reads need to be 
sorted according to the successful hits. This script takes the BWA output (BAM format)
and the raw read files, and distributes the reads into FASTA files ready for assembly.

If there are multiple results (for example, one for each read direction),
concatenate them prior to sorting.
"""

def get_filenames(directory):
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
            filenames.append(filenames_list[3])
    return filenames

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def read_sorting(bamfilename):
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
    mkdir_p(target)
    if single:
        outfile = open(
            os.path.join(target, "{}_interleaved.fasta".format(target)), 'a')
        outfile.write(">{}\n{}\n".format(ID1, Seq1))
        outfile.write(">{}\n{}\n".format(ID2, Seq2))
        outfile.close()
    else:
        outfile1 = open(os.path.join(target, "{}_1.fasta".format(target)), 'a')
        outfile1.write(">{}\n{}\n".format(ID1, Seq1))
        outfile2 = open(os.path.join(target, "{}_2.fasta".format(target)), 'a')
        outfile2.write(">{}\n{}\n".format(ID2, Seq2))
        outfile1.close()
        outfile2.close()


def write_single_seqs(target, ID1, Seq1):
    """Distributing targets from single-end sequencing"""
    mkdir_p(target)
    outfile = open(os.path.join(target, "{}_unpaired.fasta".format(target)),
                   'a')
    outfile.write(">{}\n{}\n".format(ID1, Seq1))
    outfile.close()


def distribute_reads(readfile, read_hit_dict, single=True):
    filenames = get_filenames(readfile)
    num_reads_to_write = len(read_hit_dict)
    if num_reads_to_write != 0:
        iterator1 = FastqGeneralIterator(open("".join([readfile, '/', filenames[0]])))
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
                    sys.stderr.write("[%-20s] %d%%" % ('=' * int(20 * j), 100 * j))
                    sys.stderr.flush()
            sys.stderr.write("\n")
            return

        elif len(filenames) == 2:
            iterator2 = FastqGeneralIterator(open("".join([readfile, '/', filenames[1]])))

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
                        version="distribute_reads_to_targets_bwa 1.3.2")
    parser.add_argument("-b", "--bamfile", action="store", dest="bamfile",
                        help="The location of the input bamfile",
                        required=True)
    parser.add_argument("-r", "--readfile", action="store", dest="readfile",
                        help="The location of the input readfile(s)",
                        required=True)
    argvs = parser.parse_args()
    return argvs


def main():
    argvs = parseArgvs()
    bamfilename = argvs.bamfile
    readfile = argvs.readfile
    # print(bamfilename, readfile)
    read_hit_dict = read_sorting(bamfilename)
    print(read_hit_dict)
    print("Unique reads with hits: {}".format(len(read_hit_dict)))
    distribute_reads(readfile, read_hit_dict, single=True)


if __name__ == "__main__": main()

