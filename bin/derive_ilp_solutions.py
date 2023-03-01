#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ILP file of RIBAP. This file contains several classes that are needed
for the generation of pairwise ILPs.

Usage:
ILP.py [options] <tsv_file>

Options:
    -h, --help                                      Show this neat help message and exit.
    -s SELFWEIGHT, --self SELFWIGHT                 Weight of a self-edge in the matching. [default: 0.0]
    -a ALPHA, --alpha ALPHA                         Alpha coefficient for the rearrangement/sequence balance. [default: 0.5]
    -x TOLERANCE, --fix_adj_tolerance TOLERANCE     Weight tolerance for fixing gene adjacencies. [default: 0.0]
    -max, --maximal                                 Force maximal matching.
    -ms, --matching_size                            Matching size in the F.O.
    -i, --indel                                     Use measurement of InDel events.
    -k, --keep                                      Keeps ILP files. WARNING: This takes a LOT of disk space. [default: False]

"""

__author__ = "Kevin Lamkiewicz"

from itertools import chain
import sys
import re
import pickle

import networkx as nx
from networkx.algorithms import bipartite
from networkx import connected_components

from docopt import docopt

from FFProblem import FFProblem
from ILPBuilder import ILPGenerator


def main():
  """
  Starting point for this script. All functions will be called and coordinated here.
  """
  param = docopt(__doc__)
  pickled_data, s, alpha, matching_size, maximal, fix_adj_tolerance, indel = parse_arguments(param)
  blastTable = read_blast_table(pickled_data)

  for pairwiseSpecies, similarities in blastTable.items():
    problem = FFProblem(pairwiseSpecies, similarities)
    ilpGen = ILPGenerator(pickled_data, pairwiseSpecies)
    ilpGen.generate_lp(problem, self_edge_cost=s, alpha=alpha, MATCHING_SIZE=matching_size, MAXIMAL_MATCHING=maximal, tolerance=fix_adj_tolerance, INDEL=indel)



def parse_arguments(param):
    s = float(param['--self'])
    alpha = float(param['--alpha'])
    matching_size = param['--matching_size']
    maximal = param['--maximal']
    fix_adj_tolerance = float(param['--fix_adj_tolerance'])
    indel = param['--indel']

    pickled_data=param['<tsv_file>']

    return [pickled_data, s, alpha, matching_size, maximal, fix_adj_tolerance, indel]


def read_blast_table(pickled_data):
  with open(pickled_data, 'rb') as inputStream:
    blastTable = pickle.load(inputStream)
  return blastTable

if __name__ == '__main__':
  main()