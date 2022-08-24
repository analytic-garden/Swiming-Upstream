#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
upstream_seqs.py - retreive upstream regions for othologous genes
author: Bill Thompson
license: GPL 3
copyright: 2022-08-13
"""
import sys
import argparse
import pickle
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
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

        parser = Parser(description='Retreive upstream regions for othologous genes.')

        parser.add_argument('-g', '--gene',
                            required = True,
                            type = str,
                            help = 'Target gene (required).')
        parser.add_argument('-t', '--ortholog_table',
                            required = True,
                            type = str,
                            help = 'Ortholog table (required).')
        parser.add_argument('-i', '--index_path',
                            required = True,
                            type = str,
                            help = 'Path to genome index files (required).')
        parser.add_argument('-o', '--output_file',
                            help="Path to output file. (required)",
                            required = True,
                            type = str)
        parser.add_argument('-l', '--seq_length',
                            help="Upstream sequence length (default = 1000)",
                            required = False,
                            default = 1000,
                            type = int)
            
        return parser.parse_args()

    parser = argparse.ArgumentParser(description = 'Retreive upstream regions for othologous genes.')
    args = ParseArgs(parser)

    return args

def upstream_sequence(gene, gene_features, length = 1000, min_length = 100):
    """"
    upstream_sequence - grab sequence data preceding the 5' end of a gene.

    Parameters
    ----------
    gene : str
        a gene name
    gene_features : GeneFeatures
        A GeneFeatures object for a genome.
    length : int, optional
        number of nucleotides to grab, by default 1000
    min_length : int, optional
        reject sequences shorter than this, by default 100

    Returns
    -------
    Bio.SeqRecord
        A SeqRecord containing the upstream sequence or None if region is too short.

    Requires
    --------
    The gene must exist in gene_features. 
    """
    feature = gene_features.features[gene_features.genes[gene]]
    location = feature.location
    if location.strand == 1:
        end = location.start - 1
        upstream_feature = gene_features.features[gene_features.genes[gene] - 1]
        start = max(upstream_feature.location.end + 1, location.start - length)
    else:
        start = location.end + 1
        upstream_feature = gene_features.features[gene_features.genes[gene] + 1]
        end = min(upstream_feature.location.start - 1, location.end + length)

    # overlapping genes or too short intergenic
    if end - start + 1 < min_length:
        return None

    record = SeqRecord(id = gene + '_' + gene_features.organism,
                        description = "intergenic {}:{} length {}". format(start, end, end - start + 1),
                        seq = gene_features.genome_sequence[start:(end+1)])

    return record

def main():
    args = GetArgs()
    ortholog_table_file = args.ortholog_table
    index_path = args.index_path
    gene = args.gene.upper()
    output_file = args.output_file
    seq_length = args.seq_length

    if index_path[-1] != '/':
        index_path += '/'

    ortholog_table = pd.read_csv(ortholog_table_file, delimiter = '\t')
    genome_list = list(ortholog_table['genomic_nucleotide_accession.version'])

    upstream_seqs = []
    for genome in genome_list:
        index_file = index_path + genome + '.pkl'
        with open(index_file, 'rb') as f:
            gene_features = pickle.load(f)

        rec = upstream_sequence(gene, gene_features, seq_length)
        if rec is not None:
            upstream_seqs.append(rec)

    SeqIO.write(upstream_seqs, output_file, 'fasta')

if __name__ == "__main__":
    main()
