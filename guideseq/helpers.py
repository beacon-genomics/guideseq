"""
helpers.py

A bunch of helper functions that have specific/limited uses
"""

import os
import glob

def countNumberOfLines(filepath):
    """
    Given filepath, return number of lines in that file
    """
    with open(filepath) as f:
        count = 0
        for line in f:
            count = count + 1
    return count


def countNumberOfLinesInFolder(folder_path, extension="fastq"):
    """
    Given folder path, return total number of lines in all fastq files in that folder
    """
    count = 0
    filepaths = glob.glob(os.path.join(folder_path, "*." + extension))
    for filepath in filepaths:
        count = count + countNumberOfLines(filepath)
    return count

def distanceToGene(annotation_range, gene_range):
    """
    Given the ranges of a gene and an annotation, return the distance between the two
    Based on: https://stackoverflow.com/a/16843530/1224152
    """

    # First, sort the two ranges, such that the one with the smallest first element is x
    x, y = sorted((annotation_range, gene_range))

    # Then, if x[1] is between x[0] and y[0], and x[1] != y[0], return difference
    # Else, return 0
    if x[0] <= x[1] and x[1] < y[0]:
        return y[0] - x[1]
    else:
        return 0
