#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GeneFeatures.py - a class to hold Biopython GenBank features for named genes
author: Bill Thompson
license: GPL 3
copyright: 2022-08-23
"""
from Bio import SeqIO

class GeneFeatures(object):
    def __init__(self, genome_file):
        self._genes = {}
        self._features = []
        self._index_genes(genome_file)

    @property
    def genes(self):
        return self._genes

    @property
    def features(self):
        return self._features

    @property
    def organism(self):
        return self._organism

    @property
    def genome_sequence(self):
        return  self._genome_sequence 

    def _index_genes(self, genome_file):
        """
        index_genes - create an index of features for genes in this GenBank file
                      The index consists of a dict, self.gene, with gene name as key and index of the 
                      gene feature in the feature list, self.features,

        Parameters
        ----------
        genome_file : str
            path to GenBank file.
        """
        gb = SeqIO.read(genome_file, 'genbank')
        self._genome_sequence = gb.seq

        gene_count = 0
        for feature in gb.features:
            if feature.type == 'gene':
                if 'pseudo' not in feature.qualifiers:   # we are skipping pseudo genes
                    self._features.append(feature)
                    self._genes[feature.qualifiers['gene'][0].upper()] = gene_count
                    gene_count += 1
            elif feature.type == 'source':
                self._organism = feature.qualifiers['organism'][0].replace(' ', '_')

