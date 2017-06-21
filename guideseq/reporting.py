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
import copy
import shutil
import subprocess

import yaml
import logging
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML

from gsToBed import gsToBed
from helpers import countNumberOfLines, countNumberOfLinesInFolder, distanceToGene

BEACON_LOGO_PATH = "reporting/images/BeaconGenomicsLogo.png"

def report(manifest_path):

    """
    Manifest Parsing
    """
    with open(manifest_path, 'r') as f:
        template_vars = yaml.load(f)
    print('Parsed manifest')

    output_folder = template_vars['output_folder']

    template_vars['date'] = time.strftime("%m/%d/%Y")
    template_vars['visualization_path'] = "test/data/visualization/EMX1_identifiedOfftargets_visualization.svg"
    template_vars['beacon_logo_path'] = BEACON_LOGO_PATH
    template_vars['samples_count'] = len(template_vars['samples'])

    samples = template_vars['samples']

    if 'control' in samples:
        del samples['control']

    """
    Setup
    """
    # Create temp folder
    temp_folder = os.path.join(output_folder, 'tmp')
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)
    print('Setup complete')

    """
    GLOBAL STATS
    """
    # Calculate total reads analyzed
    consolidated_path = os.path.join(template_vars['output_folder'], "consolidated")
    template_vars['total_reads'] = countNumberOfLinesInFolder(consolidated_path)

    # Calculate total number of cleavage sites found
    identified_path = os.path.join(template_vars['output_folder'], "identified")
    total_cleaved = countNumberOfLinesInFolder(identified_path)
    control_count = countNumberOfLines("test/data/identified/control_identifiedOfftargets.txt")
    template_vars['total_cleaved'] = total_cleaved - (control_count) - (template_vars['samples_count'] - 1)

    print('Calculated global stats')


    """
    CONDITIONS
    """
    # For each listed condition, determine how many samples were run in those conditions
    for i, condition in enumerate(template_vars['conditions']):
        key = list(condition.keys())[0]
        template_vars['conditions'][i][key]['count'] = 0
        for sample in template_vars['samples']:
            if template_vars['samples'][sample]['conditions'] == key:
                template_vars['conditions'][i][key]['count'] += 1

    """
    SAMPLES
    """
    samples = copy.deepcopy(template_vars['samples'])
    total_reads_command = 'samtools view {0}'
    high_quality_filter_command = 'samtools view {0} -q 50 -F 2176'
    annotation_command = 'bedtools closest -a {0} -b {1} -d > {2}'
    dnase_command = 'bedtools intersect -a {0} -b {1}'

    for sample_name, sample in samples.items():

        #Count total number of reads
        alignment_sam = os.path.join(template_vars['output_folder'], 'aligned', sample_name + '.sam')
        tmp_command = total_reads_command.format(alignment_sam)
        total_reads = subprocess.check_output(tmp_command + ' | wc -l', shell=True)
        sample['total_reads'] = int(total_reads.strip())

        # Count number of high quality mapped reads
        tmp_command = high_quality_filter_command.format(alignment_sam)
        hqmr = subprocess.check_output(tmp_command + ' | wc -l', shell=True)
        sample['hqmr'] = int(hqmr.strip())

        # Count number of on and off-target reads
        on_target_reads = 0
        off_target_reads = 0
        identified_path = os.path.join(template_vars['output_folder'], 'identified', sample_name + '_identifiedOfftargets.txt')
        identified_table = pd.read_csv(identified_path, sep='\t', index_col=False)

        on_targets = identified_table[identified_table['Mismatches'] == 0]
        off_targets = identified_table[identified_table['Mismatches'] > 0]
        cleavage_sites = pd.concat([on_targets, off_targets])

        sample['ontarget_count'] = on_targets.sum()['bi.sum.mi']
        sample['offtarget_count'] = off_targets.sum()['bi.sum.mi']

        # Make temp BED file
        bedfile_path = os.path.join(temp_folder, sample_name + '.bed')
        gsToBed(identified_path, bedfile_path)

        # Get gene annotations
        annotation_path = os.path.join(temp_folder, sample_name + '_annotation.bed')
        tmp_command = annotation_command.format(bedfile_path, template_vars['gencode_gtf'], annotation_path)
        result = subprocess.check_output(tmp_command, shell=True)

        # Read in the annotations
        annotations = pd.read_csv(annotation_path, sep='\t', index_col=False, header=None)

        # Create the cleavage annotation table
        sample['cleavage_annotations'] = []
        for i, cleavage_site in cleavage_sites.iterrows():
            # Find the corresponding annotation
            cleavage_id = cleavage_site['#BED Chromosome']
            is_cleavage = annotations[3] == cleavage_id
            is_gene = annotations[6] == 'gene'
            cleavage_annotation = annotations[is_cleavage].head(1)

            print(cleavage_annotation)
            print(type(cleavage_annotation))
            cleavage_range = cleavage_id.split(':')[1].split('-')
            annotation_range = (int(cleavage_annotation.get_value(0,1)), int(cleavage_annotation.get_value(0,2)))
            gene_distance = distanceToGene(cleavage_range, annotation_range)

            sample['cleavage_annotations'].append({
                'cleavage_id': cleavage_id,
                'location': cleavage_site['#BED Chromosome'],
                'targeting_status': 'On target' if cleavage_site['Mismatches'] == 0 else 'Off target',
                'gs_reads': 'TODO',
                'closest_gene': 'gene name',
                'distance_closest_gene': gene_distance})

            



    template_vars['samples'] = samples

    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template("reporting/templates/reportTemplate.html")

    # Generate the HTML and save it
    html_out_path = ""
    html_out = template.render(template_vars)
    # with open(html_out_path, 'wb') as f:
    #     f.write(html_out)

    # Render the HTML into PDF and save it
    pdf_out_path = "beacon_test.pdf"
    HTML(string=html_out, base_url=".").write_pdf(pdf_out_path, stylesheets=["reporting/styles/beaconReports.css"])

    """
    TEARDOWN
    """
    # Delete temp folder
    # if os.path.exists(temp_folder):
    #     shutil.rmtree(temp_folder, ignore_errors=True)

def main():
    report(sys.argv[1])

if __name__ == "__main__":
    main()
