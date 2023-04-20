#!/usr/bin/env python

'''
This .py file contains completed functions that, due to time constrains, remain unimplemented.
This file serves to document these otherwise useful functions so they may one day be implemented if 
deemed nessesary.
These functions were meant to be included in generate_hybpiper_commands.py

It includes:

A function to generate a hybpiper stats command

A function to use the stats file to generate a hybpiper recovery_heatmap command.
'''
import os, argparse


def construct_stats_command(sequence_type, targetfile_dna, targetfile_aa, namelist, seq_lengths_filename='seq_lengths', stats_filename='hybpiper_stats', run_profiler=False):
    """The function generates the hybpiper stats command by appending a template string with the proper flags according
    to the state of the arguments.

    Parameters
    ----------
    sequence_type : str
        Sequence type (gene or supercontig) to recover lengths for.
    targetfile_dna : str
        A string resembling the path to the DNA targets fasta file.
    targetfile_aa : str
        A string resembling the path to the amino-acid targets fasta file.
    namelist : str
        A text file with names of hybpiper output directories, one per line.
    seq_lengths_filename : str, optional
        File name for the sequence lengths *.tsv file. Default is <seq_lengths.tsv>.
    stats_filename : str, optional
        File name for the stats *.tsv file. Default is <hybpiper_stats>.
    run_profiler : bool, optional
        If supplied, run the subcommand using cProfile. Saves a *.csv file of results.

    Returns
    -------
    stats_cmd : str
        The generated hybpiper stats command in the form of a string.
    """
    stats_cmd = "hybpiper stats %s" % sequence_type

    if targetfile_dna:
        stats_cmd += " -t_dna %s" % targetfile_dna

    if targetfile_aa:
        stats_cmd += " -t_aa %s" % targetfile_aa

    stats_cmd += " %s" % namelist

    if seq_lengths_filename:
        stats_cmd += " --seq_lengths_filename %s" % seq_lengths_filename

    if stats_filename:
        stats_cmd += " --stats_filename %s" % stats_filename

    if run_profiler:
        stats_cmd += " --run_profiler"

    return stats_cmd, stats_filename


def construct_heatmap_command(seq_lengths_file, heatmap_filename='recovery_heatmap', figure_length=None,
                              figure_height=None, sample_text_size=None, gene_text_size=None, heatmap_filetype='png',
                              heatmap_dpi=150, run_profiler=False):
    """
    Constructs the command to create a gene recovery heatmap for a hybpiper run.

    Args:
        seq_lengths_file (str): Filename for the seq_lengths file (output of the 'hybpiper stats' command).
        heatmap_filename (str): Filename for the output heatmap, saved by default as a *.png file.
                                Defaults to "recovery_heatmap".
        figure_length (int): Length dimension (in inches) for the output heatmap file. Default is automatically
                             calculated based on the number of genes.
        figure_height (int): Height dimension (in inches) for the output heatmap file. Default is automatically
                             calculated based on the number of samples.
        sample_text_size (int): Size (in points) for the sample text labels in the output heatmap file.
                                Default is automatically calculated based on the number of samples.
        gene_text_size (int): Size (in points) for the gene text labels in the output heatmap file.
                              Default is automatically calculated based on the number of genes.
        heatmap_filetype (str): File type to save the output heatmap image as. Default is *.png.
        heatmap_dpi (int): Dot per inch (DPI) for the output heatmap image. Default is 150.
        run_profiler (bool): If True, run the subcommand using cProfile. Saves a *.csv file of results.

    Returns:
       heatmap_cmd (str): The command to create the gene recovery heatmap.
    """
    heatmap_cmd = 'hybpiper recovery_heatmap {} --heatmap_filename {} --heatmap_filetype {} --heatmap_dpi {}'.format(
        seq_lengths_file, heatmap_filename, heatmap_filetype, heatmap_dpi)

    if figure_length is not None:
        heatmap_cmd += ' --figure_length {}'.format(figure_length)
    if figure_height is not None:
        heatmap_cmd += ' --figure_height {}'.format(figure_height)
    if sample_text_size is not None:
        heatmap_cmd += ' --sample_text_size {}'.format(sample_text_size)
    if gene_text_size is not None:
        heatmap_cmd += ' --gene_text_size {}'.format(gene_text_size)
    if run_profiler:
        heatmap_cmd += ' --run_profiler'

    return heatmap_cmd


def main():
    target_format = "DNA"
    heatmap_bool = "y"

    if target_format == "DNA":
        targetfile_dna = targetfile
    elif target_format == "AA":
        targetfile_aa = targetfile
    stats_cmd, stats_filename = construct_stats_command("gene", targetfile_dna, targetfile_aa, namelist_location)
    heatmap_cmd = ""
    if heatmap_bool == "y":
        heatmap_cmd = construct_heatmap_command(stats_filename)


if __name__ == '__main__':
    main()