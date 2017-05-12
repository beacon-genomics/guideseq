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