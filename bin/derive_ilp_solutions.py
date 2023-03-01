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

    -p PROCESSES, --processes PROCESSES             Number of processes that are run in parallel to calculate
                                                    and solve ILPs [default: 1]
    -k, --keep                                      Keeps ILP files. WARNING: This takes a LOT of disk space. [default: False]

    -t TIME, --tmlim TIME                           Time (in seconds) that is allocated for each individual ILP. [default: 240]

"""

__author__ = "Kevin Lamkiewicz"

from itertools import chain
import sys
import re
import pickle
from glob import glob

from collections import defaultdict

from multiprocessing import Pool
import subprocess

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
  pickled_data, s, alpha, matching_size, maximal, fix_adj_tolerance, indel, cpus = parse_arguments(param)
  blastTable = read_blast_table(pickled_data)
  pairwiseSimple = defaultdict(list)
  
  #! TODO: Multiprocessing via pool
  for pairwiseSpecies, similarities in blastTable.items():
    problem = FFProblem(pairwiseSpecies, similarities)
    ilpGen = ILPGenerator(pickled_data, pairwiseSpecies)
    ilpGen.generate_lp(problem, self_edge_cost=s, alpha=alpha, MATCHING_SIZE=matching_size, MAXIMAL_MATCHING=maximal, tolerance=fix_adj_tolerance, INDEL=indel)
    for ilpFile in glob(f"{ilpGen.out}*ilp"):
      command = f"glpsol --lp {ilpFile} --mipgap 0.01 --pcost --cuts --memlim 16834 --tmlim 240 -o /dev/stdout"
      process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
      stdout,stderr = process.communicate()
      
      
      conditions = []
      allLine = stdout.split("\n")  
      for idx,line in enumerate(allLine):   
        if "x_A" in line:
          if len(line.split())==2:
            tmp = line + allLine[idx+1]
          else:
            tmp = line
          conditions.append(tmp)
          
      for line in conditions:
        if line.split()[3] == "1":
          pairwiseSimple[pairwiseSpecies].append(line)
          


def parse_arguments(param):
    s = float(param['--self'])
    alpha = float(param['--alpha'])
    matching_size = param['--matching_size']
    maximal = param['--maximal']
    fix_adj_tolerance = float(param['--fix_adj_tolerance'])
    indel = param['--indel']

    pickled_data=param['<tsv_file>']
    cpus = int(param['--processes'])

    return [pickled_data, s, alpha, matching_size, maximal, fix_adj_tolerance, indel, cpus]


def read_blast_table(pickled_data):
  with open(pickled_data, 'rb') as inputStream:
    blastTable = pickle.load(inputStream)
  return blastTable

if __name__ == '__main__':
  main()