#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
nice docstring
"""

__author__ = "Kevin Lamkiewicz"

import glob
import sys
import os
from collections import defaultdict
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
import shutil


dirPath = sys.argv[1]
ribapTable = sys.argv[2]
wasRenamed = sys.argv[3]

msaPath = f"{dirPath}/msa/"

id2strain = {}
aminoSequences = {}
NUMSTRAINS = 0
#for file in glob.glob(f"{dirPath}/prokka/*/*faa"):
for file in glob.glob("**/*.faa"):
  NUMSTRAINS += 1
  basename = os.path.basename(file)
  for record in SeqIO.parse(file, 'fasta'):
    if len(wasRenamed) > 1:
      basename = basename.replace(".faa", '')
    else:
      basename = basename.replace("_RENAMED.faa", '')
    geneID = record.id.split('_')[1]
    geneID = f"{basename}_{geneID}"
    id2strain[record.id] = geneID
    aminoSequences[record.id] = record.seq

group2geneIDs = {}
coreMSA = defaultdict(str)

with open(ribapTable, 'r') as inputStream:
  for line in inputStream:
    if line.startswith("Cluster"): continue
    array = line.rstrip().split('\t')
    group2geneIDs[array[0]] = array[3:]
    

coreGroups = []

for group, genes in group2geneIDs.items():
  genes = [x for x in genes if x != 'NA']
  sequences = [aminoSequences[x] for x in genes]
  if len(sequences) <= 2:
    continue
  
  fastaRecords = zip(genes, sequences)
  with open(f"{msaPath}/{group}.faa", 'w') as outputStream:
    for geneID, seq in fastaRecords:
      outputStream.write(f">{id2strain[geneID]}\n{seq}\n")

      if len(sequences) == NUMSTRAINS:
        coreGroups.append(group)
        #strain = '_'.join(id2strain[geneID].split("_")[:-1])
        #coreMSA[strain] += seq

# coreMSA = defaultdict(str)
# 
#for file in glob.glob(f"{msaPath}/*.faa"):
  # mafftBin = shutil.which('mafft')
  # cmd = MafftCommandline(mafftBin, input=file)
  # stdout, stderr = cmd()
  # alnOut = file.replace('.faa', '')
  # 
  # with open(f'{alnOut}_mafft.aln', 'w') as outputStream:
    # outputStream.write(stdout)
  # if stdout.count(">") == NUMSTRAINS:
    # header = ""
    # sequence = ""
    # for line in stdout.split("\n"):
      # if line.startswith(">"):
        # if header:
          # coreMSA[header] += sequence
        # header = "_".join(line.split('_')[:-1])
        # sequence = ""
        # continue
      # sequence += line
    # coreMSA[header] += sequence
# 
# with open(f"{msaPath}/coreGenome_mafft.aln", 'w') as outputStream:
  # for header, sequence in coreMSA.items():
    # outputStream.write(f"{header}\n{sequence}\n")
    