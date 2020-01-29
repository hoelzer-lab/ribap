#!/usr/bin/env python3

import sys
import os
import glob

dirPath = sys.argv[1]
NUMSTRAINS = int(sys.argv[2])

coreAlignment = {}

for file in glob.glob(f"{dirPath}/*aln"):
  with open(file, 'r') as inputStream:
    
    strainsInFile = sum([1 for line in inputStream if line.startswith(">")])
    if not strainsInFile == NUMSTRAINS:
      continue
    inputStream.seek(0)

    strainID = ""
    sequence = ""

    for line in inputStream:
      if line.startswith(">"):
        if strainID:
          if strainID in coreAlignment.keys():
          #if strainID:
            coreAlignment[strainID] +=  "XXXXX"+sequence          
          else:
            coreAlignment[strainID] = sequence          
        strainID = "".join((line.lstrip('>').split('_')[:-1]))
        sequence = ""
      else:
        sequence += line.rstrip()

    if strainID in coreAlignment.keys():
      coreAlignment[strainID] += "XXXXX"+sequence          
    else:
      coreAlignment[strainID] = sequence          


with open(f"{dirPath}/coreGenome_mafft.aln", 'w') as outputStream:
  for strainID, sequence in coreAlignment.items():
    outputStream.write(f">{strainID}\n{sequence}\n")

