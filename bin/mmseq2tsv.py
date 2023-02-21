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

# This can't happen anymore. If it does, we screwed up!
if len(sys.argv) != 4:
    print("Not enough arguments.")
    print("Are you calling this script yourself, even though you are not supposed to?")
    print("Well, nothing will happen, thus RIBAP is exiting.")
    sys.exit(1)

blast = sys.argv[1]
strainID = sys.argv[2]
outputPath = sys.argv[3]


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
#combs = itertools.combinations(strains.keys(), 2)
#for i in combs:
#    strains = sorted(i)
#    key = f"{strains[0]}:{strains[1]}"
#    stream = open(f"{outputPath}/{strains[0]}-vs-{strains[1]}.tsv", 'w')
#    outputStreams[key] = stream

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
        #orientation = -1 if int(currentEntry[8]) > int(currentEntry[9]) else 1
        # with open(f"{outputPath}/{queryStrain}-vs-{targetStrain}.tsv", 'a') as outputStream:
        #     outputStream.write(f"{queryGene}\t{targetGene}\t{orientation}\t{seqSim}\n")

        #key = f"{queryID.split('_')[0]}:{targetID.split('_')[0]}"
        #key = f"{'_'.join(queryID)}:{'_'.join(targetID)}"
        #stream = outputStreams[key]
        #stream.write(f"{int(queryID.split('_')[1])}\t{int(targetID.split('_')[1])}\t{orientation}\t{seqSim}\n")
        
        key = (queryStrain, targetStrain)

        if key not in blastTable:
            blastTable[key] = [(queryGene, targetGene, seqSim, orientation)]
        else:
            blastTable[key].append((queryGene, targetGene, seqSim, orientation))
# exit(0)
#for stream in outputStreams.values():
#    stream.close()

with open("mmseqs_compressed.pkl", 'wb') as f:
    pickle.dump(blastTable, f)

#sys.exit(0)
# for each strain pairwise comparison, write a file with all gene entries
#pairwiseComparisons = []

# for key, info in blastTable.items():
#     print(key, info)
#     seqSim = info[0]
#     orientation = info[1]
#     # lots of parsing from our hashes, nothing super special
#     queryID = key.split(':')[0].split("_")
#     targetID = key.split(':')[1].split("_")
#     #print(queryID)
#     queryStrain = queryID[0]
#     queryGene = int(queryID[1])
#     targetStrain = targetID[0]
#     targetGene = int(targetID[1])
#     bitscore = int(seqSim)

    #with open(f"{outputPath}/{queryStrain}-vs-{targetStrain}.tsv", 'a') as outputStream:
    #    outputStream.write(f"{queryGene}\t{targetGene}\t{orientation}\t{seqSim}\n")



# for strainA in strains.keys():
#     for strainB in strains.keys():
# #        # again, self-hits are ignored
#         if strainA == strainB:
#             continue
# #
# #        # check whether we already have that file
#         if strainA < strainB:
#             first = strainA
#             second = strainB
#         else:
#             first = strainB
#             second = strainA
            
#         blastKey = f"{first}:{second}"
#         blastKeyRev = f"{second}:{first}"
#         if blastKey in pairwiseComparisons or blastKeyRev in pairwiseComparisons:
#             continue
# #        
#         pairwiseComparisons.append(blastKey)
#         pairwiseComparisons.append(blastKeyRev)

# print(pairwiseComparisons)
# #        
#         with open(f"{outputPath}/{first}-vs-{second}.tsv", 'w') as outputStream:
#             for key, info in blastTable.items():
#                 seqSim = info[0]
#                 orientation = info[1]
#                 # lots of parsing from our hashes, nothing super special
#                 queryID = key.split(':')[0]
#                 targetID = key.split(':')[1]
#                 queryStrain = queryID.split('_')[0]
#                 queryGene = int(queryID.split('_')[1])
#                 targetStrain = targetID.split('_')[0]
#                 targetGene = int(targetID.split('_')[1])
#                 # if the blast-key is the same as the current strains we are looking at, then we have to write
#                 # otherwise simply continue with the next blast hit
#                 if strainA == queryStrain and strainB == targetStrain:
#                     outputStream.write(f"{queryGene}\t{targetGene}\t{orientation}\t{seqSim}\n")

sys.exit(0)