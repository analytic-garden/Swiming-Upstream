#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
index_genomes.py - Build an index of genes and save it to file
author: Bill Thompson
license: GPL 3
copyright: 2022-08-23
"""
import sys
import os
import argparse
import pickle
from pathlib import Path
from GeneFeatures import GeneFeatures

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

        parser = Parser(description='Build an index of genes and save it to file.')

        parser.add_argument('-p', '--genome_path',
                            required = True,
                            type = str,
                            help = 'Path to genome files (required).')
        parser.add_argument('-o', '--output_path',
                            help="Path to output files. (default = geneom_path)",
                            required = False,
                            default = None,
                            type = str)
            
        return parser.parse_args()

    parser = argparse.ArgumentParser(description = 'Build an index of genes and save it to file.')
    args = ParseArgs(parser)

    return args

def main():
    args = GetArgs()
    genome_path = args.genome_path
    output_path = args.output_path

    if genome_path[-1] != '/':
        genome_path += '/'
    if output_path is None:
        output_path = genome_path
    elif output_path[-1] != '/':
        output_path += '/'

    for p in Path(genome_path).glob('*.gb'):
        features = GeneFeatures(genome_path + p.name)
        out_file = output_path + os.path.basename(os.path.splitext(p.name)[0]) + '.pkl'
        with open(out_file, 'wb') as f:
            pickle.dump(features, f)

if __name__ == "__main__":
    main()
