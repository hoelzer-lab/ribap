#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file a small script in the RIBAP pipeline.
It is not supposed to be called directly, but from the RIBAP workflow itself.

The roary blast output is parsed into several .tsv tables that serve as input
for our ILP generator.
"""

__author__ = "Kevin Lamkiewicz"

import sys
import itertools
import pickle


def chunks(data, size):
    iterator = iter(data)
    for i in range(0, len(data), size):
        yield {k:data[k] for k in itertools.islice(iterator, size)}


# This can't happen anymore. If it does, we screwed up!
if len(sys.argv) != 5:
    print("Not enough arguments.")
    print("Are you calling this script yourself, even though you are not supposed to?")
    print("Well, nothing will happen, thus RIBAP is exiting.")
    sys.exit(1)

blast = sys.argv[1]
strainID = sys.argv[2]
outputPath = sys.argv[3]
cores = int(sys.argv[4])


# Reading strain IDs and save it to a hash.
strains = {}
with open(strainID, 'r') as inputStream:
    while 1:
        currentLine = inputStream.readline()
        if not currentLine: # EOF
            break
        strain = currentLine.rstrip("\n").split(",")
        strains[strain[0]] = strain[1]

outputStreams = {}

# Read blast table and save pairwise hits in hash.
blastTable = {}
with open(blast, 'r') as inputStream:
    while 1:
        currentLine = inputStream.readline()
        if not currentLine: # EOF
            break
        currentEntry = currentLine.rstrip("\n").split()

        queryID = currentEntry[0]
        targetID = currentEntry[1]
        if queryID > targetID:
            targetID = currentEntry[0]
            queryID = currentEntry[1]
        
        queryID = queryID.split('_')
        queryStrain = queryID[0]
        queryGene = int(queryID[1])
        targetID = targetID.split('_')
        targetStrain = targetID[0]
        targetGene = int(targetID[1])

        # if two genes from the same strain have a blast hit
        # we skip that one.
        if queryStrain == targetStrain:
            continue

        seqSim = currentEntry[-1]

        # orientation is needed to encode whether a gene is reversed in a strain
        if int(currentEntry[6]) > int(currentEntry[7]):
            orientation = -1
        else:
            orientation = 1    
        
        key = (queryStrain, targetStrain)

        if key not in blastTable:
            blastTable[key] = [(queryGene, targetGene, seqSim, orientation)]
        else:
            blastTable[key].append((queryGene, targetGene, seqSim, orientation))

if not blastTable:
    sys.exit('blastTable dictionary is empty, please check mmseqs2 output TSV is not empty and there is at least one queryStrain != targetStrain.')

print(blastTable)
print(len(blastTable), cores)
print(len(blastTable) / cores)
chunksize = int(len(blastTable) / cores) if int(len(blastTable) / cores) >= 1 else 1


for idx, item in enumerate(chunks(blastTable, chunksize)):
    with open(f"mmseqs_compressed_chunk{idx}.pkl", 'wb') as f:
        pickle.dump(item, f)


sys.exit(0)