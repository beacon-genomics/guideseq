"""
report.py

Take the results of a GUIDE-Seq run and produce a results report.

The following are required inputs:
- Dictionary of manifest variables
- Path to folder of GUIDE-Seq outputs
- Path to the refFlat(?)
"""

import os
import sys
import glob
import time 

import yaml
import logging
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML

BEACON_LOGO_PATH = "reporting/images/BeaconGenomicsLogo.png"

# logger = logging.getLogger('root')

def countNumberOfLines(filepath):
    """
    Given filepath, return number of lines in that file
    """
    with open(filepath) as f:
        count = 0
        for line in f:
            count = count + 1
    return count


def countNumberOfLinesInFolder(folder_path):
    """
    Given folder path, return total number of lines in all files in that folder
    """
    count = 0
    filepaths = glob.glob(os.path.join(folder_path, "*.fastq"))
    for filepath in filepaths:
        count = count + countNumberOfLines(filepath)
    return count


def report(manifest_path):
    with open(manifest_path, 'r') as f:
        template_vars = yaml.load(f)

    output_folder = template_vars['output_folder']

    template_vars['date'] = time.strftime("%m/%d/%Y")
    template_vars['visualization_path'] = "test/reporting/visualization/EMX1_identifiedOfftargets_visualization.svg"
    template_vars['beacon_logo_path'] = BEACON_LOGO_PATH
    template_vars['samples_count'] = len(template_vars['samples'])

    # For each sample, determine everything...
    for sample in template_vars['samples']:
        if sample == "control":
            pass
        

    # Calculate total reads analyzed
    consolidated_path = os.path.join(template_vars['output_folder'], "consolidated")
    template_vars['total_reads'] = countNumberOfLinesInFolder(consolidated_path)

    # Calculate total number of cleavage sites found
    identified_path = os.path.join(template_vars['output_folder'], "identified")
    total_cleaved = countNumberOfLinesInFolder(identified_path)
    control_count = countNumberOfLines("test/reporting/identified/control_identifiedOfftargets.txt")
    template_vars['total_cleaved'] = total_cleaved - (control_count) - (template_vars['samples_count'] - 1)

    # Calculate number of on-target sites
    

    # Calculate total number of off-target sites

    # Determine 

    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template("reporting/template/reportTemplate.html")

    # Generate the HTML and save it
    html_out_path = ""
    html_out = template.render(template_vars)
    # with open(html_out_path, 'wb') as f:
    #     f.write(html_out)

    # Render the HTML into PDF and save it
    pdf_out_path = "beacon_test.pdf"
    HTML(string=html_out, base_url=".").write_pdf(pdf_out_path, stylesheets=["reporting/style/beaconReports.css"])

def main():
    report(sys.argv[1])

if __name__ == "__main__":
    main()
