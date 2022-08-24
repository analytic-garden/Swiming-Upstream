#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fetch_genomes_from ncbi.py - retreive genome files from NCBI
author: Bill Thompson
license: GPL 3
copyright: 2022-08-11

Genome IDs are contained in an ortholog table downloaded from NCBI gene.
"""
import sys
import argparse
import pandas as pd
from Bio import Entrez
from Bio import SeqIO

def GetArgs():
    """
    GetArgs - return arguments from command line. 
    Use -h to see all arguments and their defaults.

    Returns
    args - Parameter values.
    TYPE A Parser object.
    """
    def ParseArgs(parser):
        class Parser(argparse.ArgumentParser):
            def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)

        parser = Parser(description='Retreive GenBank genome files from NCBI.')

        parser.add_argument('-t', '--ortholog_table',
                            required = True,
                            type = str,
                            help = 'Ortholog table from NCBI gene (required).')
        parser.add_argument('-g', '--genome_path',
                            help="Path to output files. (required)",
                            required = True,
                            type = str)
            
        return parser.parse_args()

    parser = argparse.ArgumentParser(description = 'Retreive GenBank genome files from NCBI.')
    args = ParseArgs(parser)

    return args

def main():
    args = GetArgs()
    genome_path = args.genome_path
    ortholog_table_file = args.ortholog_table

    if genome_path[-1] != '/':
        genome_path += '/'

    # read the table and get a list of genomes
    ortholog_table = pd.read_csv(ortholog_table_file, delimiter = '\t')
    genome_list = ortholog_table['genomic_nucleotide_accession.version'].to_list()

    Entrez.email = 'analyticgarden@gmail.com'
    Entrez.tool = 'BioPython'

    # get the file from NCBI and write it in GenBank format
    for genome_id in genome_list:
        handle = Entrez.efetch(db = 'nuccore', id = genome_id, rettype = 'gbwithparts', retmode = 'text')
        record = SeqIO.read(handle, format = 'genbank')
        filename = genome_path + genome_id + '.gb'
        print('Saving', filename)
        SeqIO.write(record, filename, 'genbank')

if __name__ == "__main__":
    main()
