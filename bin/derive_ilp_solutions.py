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

    -k, --keep                                      Keeps ILP files. WARNING: This takes a LOT of disk space.

    -t TIME, --tmlim TIME                           Time (in seconds) that is allocated for each individual ILP. [default: 240]

"""

__author__ = "Kevin Lamkiewicz"

from itertools import chain
import sys
import os
import re
import pickle
from glob import glob

from collections import defaultdict

from multiprocessing import Pool, Manager
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
  global pickled_data
  global s
  global alpha
  global matching_size
  global maximal
  global fix_adj_tolerance
  global indel
  global tmlim
  global shes_a_keeper
  global pairwiseSimple

  param = docopt(__doc__)
  pickled_data, s, alpha, matching_size, maximal, fix_adj_tolerance, indel, tmlim, shes_a_keeper = parse_arguments(param)
  blastTable = read_blast_table(pickled_data)
  

  manager = Manager()
  pairwiseSimple = manager.dict()

  for pairwiseSpecies, similarities in blastTable.items():
    pool_workload(pairwiseSpecies, similarities)


def parse_arguments(param):
    s = float(param['--self'])
    alpha = float(param['--alpha'])
    matching_size = param['--matching_size']
    maximal = param['--maximal']
    fix_adj_tolerance = float(param['--fix_adj_tolerance'])
    indel = param['--indel']

    pickled_data=param['<tsv_file>']
    tmlim = int(param['--tmlim'])
    shes_a_keeper = param['--keep']

    return [pickled_data, s, alpha, matching_size, maximal, fix_adj_tolerance, indel, tmlim, shes_a_keeper]


def read_blast_table(pickled_data):
  """

  """

  with open(pickled_data, 'rb') as inputStream:
    blastTable = pickle.load(inputStream)
  return blastTable


def pool_workload(pairwiseSpecies, similarities):
  """

  """
  #global pairwiseSimple
  #pairwiseSpecies, similarities = x
  # pairwiseSimple[pairwiseSpecies] = []
  problem = FFProblem(pairwiseSpecies, similarities)
  ilpGen = ILPGenerator(pickled_data, pairwiseSpecies)
  ilpGen.generate_lp(problem, self_edge_cost=s, alpha=alpha, MATCHING_SIZE=matching_size, MAXIMAL_MATCHING=maximal, tolerance=fix_adj_tolerance, INDEL=indel)
  conditions = []

  for ilpFile in glob(f"{ilpGen.out}*ilp"):
    command = f"glpsol --lp {ilpFile} --mipgap 0.01 --pcost --cuts --memlim 16834 --tmlim {tmlim} -o /dev/stdout"
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
    stdout,stderr = process.communicate()
    
    allLine = stdout.split("\n")  
    for idx,line in enumerate(allLine):   
      if "x_A" in line:
        if len(line.split())==2:
          tmp = line + allLine[idx+1]
        else:
          tmp = line
        conditions.append(tmp)
    if not shes_a_keeper:
      os.remove(ilpFile)

  with open(f"{pairwiseSpecies[0]}-vs-{pairwiseSpecies[1]}.ilp.simple", 'w') as outputStream:
    for line in conditions:
      if line.split()[3] == "1":  
        outputStream.write("".join(line)+"\n")   
      #pairwiseSimple[pairwiseSpecies] = pairwiseSimple[pairwiseSpecies] + [line]


# def write_simple_solutions(pairwiseSimple):
#   """

#   """
#   for pairwiseSpecies, simpleRelations in pairwiseSimple.items():
#     with open(f"{pairwiseSpecies[0]}-vs-{pairwiseSpecies[1]}.ilp.simple", 'w') as outputStream:
#       for simpleRelation in simpleRelations:
#         outputStream.write("".join(simpleRelation)+"\n")



if __name__ == '__main__':
  main()