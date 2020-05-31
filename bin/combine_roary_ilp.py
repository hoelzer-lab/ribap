#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file a small script in the RIBAP pipeline.
It is not supposed to be called directly, but from the RIBAP workflow itself.

Roary and ILP outputs are combined to create a convenient overview of the
core and pangenome.
"""

__author__ = "Kevin Lamkiewicz"

import sys
import os
import glob
import re
from collections import defaultdict
from collections import Counter

def adjust_length(identifier):
    """
    """
    name = identifier.split('_')[0]
    count = identifier.split('_')[1]
    fillZero = '0' * ( 5 - len(count))
    return(f"{name}_{fillZero}{count}")

def read_strains(strainFile):
    """
    """
    # ToDo: why local variables?
    s2i = {}
    i2s = {}
    with open(strainFile, 'r') as inputStream:
        for currentLine in inputStream:
            currentLine = currentLine.rstrip('\n')
            s2i[currentLine.split(',')[1]] = currentLine.split(',')[0]
            i2s[currentLine.split(',')[0]] = currentLine.split(',')[1]
    return s2i, i2s

def read_roary_table(roaryFile):
    """
    """
    unsplitParalogs = []
    c = 0
    with open(roaryFile, 'r') as inputStream:
        for line in inputStream:
        #while 1:
            currentLine = line.rstrip('\n')
         #   if not currentLine:
         #       break
        
            # split for "," because Prokka puts ',' into the annotation
            # and because roary puts '"' into the output table
            # which is a pain in the ass.
            array = currentLine.split('","')
            formattedArray = []
            for entry in array:
                formattedArray.append(f'"{entry}"') # ToDo: why reintroduce the double quotes?
        
            # header line, just some strain parsing
            if currentLine.startswith('"Gene"'):
                for column in POS:
                    strain = formattedArray[column].replace('"','').strip()
                    strainsRoary.append(strain)
            else:
                c += 1
                clusterSizes[f"cluster{c}"] = 0
                for column in POS:
                    geneID = formattedArray[column].replace('"', '').strip()
                    # if a strain has a gene for the corresponding roary cluster
                    if geneID:
                        # roary puts tab as separator for paralog genes.
                        if '\t' in geneID:
                            geneID, paralogs = geneID.split('\t')[0], geneID.split('\t')[1:]
                            unsplitParalogs.extend(paralogs) 
                        clusterSizes[f"cluster{c}"] += 1
                        # sanity check. This should never happen
                        if geneID in genes.keys():
                            sys.exit(1)
                        # for each gene from each strain, save the cluster it belongs to
                        genes[geneID] = f"cluster{c}"
                        # and for each cluster save list of genes
                        cluster2gene[f"cluster{c}"].append(geneID)

        #ToDo: paralog splitting may mess up cluster
        for paralog in unsplitParalogs:
            c += 1
            clusterSizes[f"cluster{c}"] = 0
            genes[paralog] = f"cluster{c}"
            clusterSizes[f"cluster{c}"] += 1
            cluster2gene[f"cluster{c}"].append(paralog)
        #print(unsplitParalogs)

def read_pairwise_ILPs(ilpPath):
    """
    """
    for ilp in glob.glob(f"{ilpPath}/*simple"):

        basename = os.path.basename(ilp).split('.')[0] 
        strain1 = basename.split('-vs-')[0]
        strain2 = basename.split('-vs-')[1]
        comparison1 = f"{strain1}:{strain2}"
        comparison2 = f"{strain2}:{strain1}"
        ilps[comparison1] = {}
        ilps[comparison2] = {}
        geneRelations = {}
        with open(ilp, 'r') as inputStream:
            for line in inputStream:
                currentLine = line.rstrip("\n")
                # store the ILP edge (x_A1_B1) in geneRelation
                geneRelation = currentLine.split()[1]
                # extract gene numbers from relation
                geneA = geneRelation.split('_')[1]
                geneB = geneRelation.split('_')[2]
                # for each gene, store the relation
                #ToDo: is geneA already key of geneRelations?
                #if geneA in geneRelations:
                #    print(f"No way! {geneA} is already here!")
                geneRelations[geneA] = geneRelation
            

            # for each gene in Strain 1, get the relation of both
            # gene extremities and check for consistency
            for headInA, relationHeads in geneRelations.items():
                # only check for one extremity per gene
                # because the ILP enforces the consistency
                if not headInA.endswith('h'):
                    continue

                tailInA = headInA.replace('h','t')
                if tailInA == geneRelations[headInA].split('_')[2]:
                    continue
                relationTails = geneRelations[tailInA]

                # weird parsing stuff;  x_A456h_B458h --> DCKEFHG_456 and TYQFDEEW_458
                # map ILP syntax back to Prokka gene IDs
                # A2h_B4t
                headGene1 = relationHeads.split('_')[1].replace('h','').replace('A', f"{strain1}_")
                headGene2 = relationHeads.split('_')[2].replace('t','').replace('h','').replace('B', f"{strain2}_")
                #headGene2 = relationHeads.split('_')[2].split('h')[0].replace('B', f"{strain2}_")
                tailGene1 = relationTails.split('_')[1].replace('t','').replace('A', f"{strain1}_")
                tailGene2 = relationTails.split('_')[2].replace('h','').replace('t','').replace('B', f"{strain2}_")
                #tailGene2 = relationTails.split('_')[2].split('t')[0].replace('B', f"{strain2}_")

                # call adjust_length for be consistent with Prokka gene ID
                # DCKEFHG_456 and TYQFDEEW_458 --> DCKEFHG_00456 and TYQFDEEW_00458
                headGene1 = adjust_length(headGene1)
                headGene2 = adjust_length(headGene2)
                tailGene1 = adjust_length(tailGene1)
                tailGene2 = adjust_length(tailGene2)
                # check for consistency in ILPs
                if headGene1 == tailGene1 and headGene2 == tailGene2:
                    # if everything is fine, store the homology in the ilp dict
                    ilps[comparison1][headGene1] = headGene2
                    ilps[comparison2][headGene2] = headGene1
                else:
                    #FixMe: Write a custom error that gets thrown here.
                    print(f"Holy Guacamolee! We lost a gene! {headGene1}, {tailGene1}, {headGene2}, {tailGene2}")

def read_prokka_annotions(path):
    """
    """
    regexAnno = re.compile(r'product=(.*)')
    regexName = re.compile(r'Name=(\w*);?')
    regexID = re.compile(r'ID=(\w{8}_\d{5})')
    #annoFiles = [x for x in glob.glob(f"{path}/../prokka/*/*.gff")]
    annoFiles = [x for x in glob.glob("**/*.gff")]

    geneAnnotations = {}
    geneNames = {}

    for file in annoFiles:
        with open(file, 'r') as inputStream:
            for line in inputStream:
                if line.startswith('##'):
                    continue
                try:
                    geneID = regexID.findall(line)[0]
                    annotation = regexAnno.findall(line)[0]
                    name = regexName.findall(line)
                    
                    if not name:
                        name = '--'
                    else:
                        name = name[0]                
                except IndexError:
                    continue # these are the sequence information lines in the prokka annotation
                
                geneAnnotations[geneID] = annotation
                geneNames[geneID] = name
    
    return geneAnnotations, geneNames

def create_ribap_groups():
    """
    """
    groupCounter = 0
    tmpGcounter = 0

    for cluster, geneList in cluster2gene.items():
  
        if any([cluster in ribapGroup for ribapName, ribapGroup in assignedGroups.items()]):
            continue

        clusters = set([cluster])

        for gene in geneList:
            strain = gene.split('_')[0]

            for comparison, results_h in ilps.items():
                strain1 = comparison.split(':')[0]
            # if the ILP comparison is with our current strain
                if strain1 == strain:
                # if in the current pairwise comparison
                # our gene of interest has a hit with the second strain
                    if gene in results_h:
                    # then get the homologGene from the second strain (based on ILP)
                        homologGene = results_h[gene]    
                        if homologGene in genes:
                        # if the ILP-based homolog gene has a assigned roary cluster
                        # merge the clusters of the current gene and homolog gene into
                        # one RIBAP group
                        # genes == gene2cluster
                            clusters.add(genes[homologGene])

            clusterDone = False
            for cl in clusters:
                if clusterDone: 
                    break
                for ribapName, ribapGroup in  assignedGroups.items():
                    if cl in ribapGroup:
                        assignedGroups[ribapName] = ribapGroup.union(clusters)
                        clusterDone = True
                        break


            if not clusterDone:
                groupCounter += 1
                assignedGroups[f'group{groupCounter}'] = set(clusters)

def merge_paralogs_to_subgroup(groupID, strain2paralogs, geneHits, subgroupCounter, subgroups={}):
    """
    """
    paralogForMainGroup = {}
    for strain, paralogs in strain2paralogs.items():
        if not paralogs:
            continue
        paralogWithMaxScore = None
        paralogFound = False

        ilpScore = 0
        ilpScores = []

        roaryScore = 0
        roaryScores = []

        ## 1) select gene based on ILP support for strain, paralogs in strain2paralogs.items():
        # which of our paralog genes has the most ILP hits to other strains
        for gene in paralogs:
            ilpScores.append(len(geneHits[gene]))
            if len(geneHits[gene]) > ilpScore:
                paralogWithMaxScore = gene
                ilpScore = len(geneHits[gene])

        ## 2) if all dudes have the same ILP support, select the dude with the best Roary support
        if len(set(ilpScores)) == 1:
            for gene in paralogs:
                # check how many genes are inside the cluster of the cluster of our current paralog gene
                roaryScores.append(len(cluster2gene[genes[gene]]))
                if len(cluster2gene[genes[gene]]) > roaryScore:
                    paralogWithMaxScore = gene
                    roaryScore = len(cluster2gene[genes[gene]])

        ## 3) if all dudes have the same ILP support and Roary score, just select the first occuring dude
        

        if len(set(ilpScores)) == 1 and len(set(roaryScores)) == 1:
            #FixMe: make a better choice based on GeneNames and/or Annotation
            currentGroup = assignedGroups[groupID]
            
            annotations = [geneAnnotations[gene] for gene in currentGroup]
            mostCommonAnnotation = max(annotations, key=annotations.count)

            for candidateParalog in paralogs:
                candidateAnno = geneAnnotations[candidateParalog]
                if candidateAnno == mostCommonAnnotation:
                    paralogFound = True
                    paralogWithMaxScore = candidateParalog
                    break
        else:
            paralogFound = True
        
        if not paralogFound:
            paralogWithMaxScore = paralogs[0]
        paralogForMainGroup[strain] = paralogWithMaxScore
        try:
            paralogs.remove(paralogWithMaxScore)
        except ValueError:
            #FixMe: provide a better error handling here
            print(groupID, subgroupCounter, paralogs, paralogWithMaxScore)
            print(strain2paralogs)
            exit(0)
    
    subgroupGenes = set([y for x in strain2paralogs.values() for y in x])
    subgroups.update({f"{groupID}.{subgroupCounter}" : subgroupGenes})

    for subID, paralogsInSub in subgroups.items():            
        if subID == f"{groupID}.{subgroupCounter}":
            continue
        for usedGene in subgroupGenes:
            if usedGene in paralogsInSub:
                paralogsInSub.remove(usedGene)

    if not all([len(x) <= 1 for x in strain2paralogs.values()]):
        subgroups.update(merge_paralogs_to_subgroup(groupID, strain2paralogs, geneHits, subgroupCounter+1, subgroups))
    
    try:
        for removableGene in subgroupGenes:
            assignedGroups[groupID].remove(removableGene)    
    except ValueError:
        exit(1)
    return subgroups

def fill_group_with_genes():
    """
    """
    subgroups = {}
    for group, clusterArray in assignedGroups.items():
        genesOfGroup = set([])
        
        genesOfGroup = set([x for cluster in clusterArray for x in cluster2gene[cluster]])
        assignedGroups[group] = genesOfGroup
        strainsInGroup = [x.split('_')[0] for x in genesOfGroup]

        # trivial line of code...
        strain2paralogs = { strain : [genes for genes in genesOfGroup if genes.startswith(strain)] for strain in strainsInGroup if strainsInGroup.count(strain) > 1 }


        if strain2paralogs:
            ilpSupport = {}
            

            for strain, paralogs in strain2paralogs.items():
                for paralogGene in paralogs:
                    ilpSupport[paralogGene] = set()
                    for remainingGene in genesOfGroup:
                        if remainingGene.startswith(strain):
                            continue
                        ilpComparison = f"{strain}:{remainingGene.split('_')[0]}"
                        # ilps[ilpComparison]: { strainA_geneID : strainB_geneID }
                        if paralogGene in ilps[ilpComparison]:
                            if ilps[ilpComparison][paralogGene] in genesOfGroup:
                                ilpSupport[paralogGene].add(ilps[ilpComparison][paralogGene])
            
            subgroups.update(merge_paralogs_to_subgroup(group, strain2paralogs, ilpSupport, 1))
                
    assignedGroups.update(subgroups)

def remove_duplicated_groups(assignedGroups):
    groups = sorted(assignedGroups, key=lambda x: len(assignedGroups[x]))
    groupsToDelete = set()

    for idx, groupID in enumerate(groups):
        for largerID in groups[idx+1:]:
            if assignedGroups[groupID].intersection(assignedGroups[largerID]):
                assignedGroups[largerID] = assignedGroups[groupID].union(assignedGroups[largerID])
                groupsToDelete.add(groupID)
    assignedGroups = {groupID : content for groupID, content in assignedGroups.items() if groupID not in groupsToDelete}
    return assignedGroups

def write_output(outputFile, ident):
    sortedStrainIDs = []
 #   annoTable = os.path.dirname(outputFile) + "/ribap_individual_annotation.csv"
    annoTable = "ribap_individual_annotation" + ident + ".csv"

    for roaryInputStrain in strainsRoary:
        sortedStrainIDs.append(strain2id[roaryInputStrain])

    with open(annoTable, 'w') as outputStreamAnno:
        with open(outputFile, 'w') as outputStream:
            # write the header line of the table
            outputStream.write("Cluster_ID\t"+"Annotation\t"+"Gene_Name\t"+'\t'.join([id2strain[x].replace("_RENAMED","") for x in sortedStrainIDs])+"\n")
            outputStreamAnno.write("Cluster_ID\t"+'\t'.join([id2strain[x].replace("_RENAMED","") for x in sortedStrainIDs])+"\n")
            
            for group, geneArray in assignedGroups.items():
                if not geneArray:
                    continue
                annotations = [geneAnnotations[x] for x in geneArray]
                names = [geneNames[x] for x in geneArray]
                
                majorityAnnotation = Counter(annotations).most_common(1)[0][0]
                majorityName = Counter(names).most_common(1)[0][0]                
                

                tmpHash = { strainID : "NA" for strainID in sortedStrainIDs }
                for strainID in sortedStrainIDs:
                    for geneID in geneArray:
                        if geneID.startswith(strainID):
                            tmpHash[strainID] = geneID
                outputRow = list(tmpHash.values())
                
                annoOutput = [geneAnnotations[x] if x in geneAnnotations else "NA" for x in outputRow]
                nameOutput = [geneNames[x] if x in geneNames else '--' for x in outputRow]
                geneDescription = [f"{name} // {desc}" for name, desc in zip(nameOutput,annoOutput)]
                
                outputStream.write(group+"\t"+majorityAnnotation+"\t"+majorityName+"\t"+"\t".join(outputRow)+"\n")
                outputStreamAnno.write(group+"\t"+"\t".join(geneDescription)+"\n")

##########################################

# structures we need during the script
genes = {}
cluster2gene = defaultdict(list)
clusterSizes = {}
strainsRoary = []
assignedGroups = {}
##########################################

# read strain IDs from Prokka
strain2id, id2strain = read_strains(sys.argv[1])
NUMSTRAINS = len(strain2id.keys())
POS = list(range(14,14+NUMSTRAINS))

##########################################

# read and parse the roary table
# function for roary parsing
read_roary_table(sys.argv[2])

##########################################

# reading pairwise ILPs
ilps = {}
read_pairwise_ILPs(sys.argv[3])

##########################################

# read prokka gene annotations
geneAnnotations, geneNames = read_prokka_annotions(sys.argv[3])

##########################################

# create ribap groups based on roary and ILPs
create_ribap_groups()
assignedGroups = remove_duplicated_groups(assignedGroups)

# assign genes to the roary clusters of each ribap group
fill_group_with_genes()
assignedGroups = dict(sorted(list(assignedGroups.items()), key=lambda x: float(x[0].replace("group",''))))

write_output(sys.argv[4], sys.argv[5])

##########################################
coreGeneCounter = 0
for group, strains in assignedGroups.items():
    if len(strains) == NUMSTRAINS:
        coreGeneCounter += 1

print(f"We summarized {len(clusterSizes.keys())} Roary clusters into {len(assignedGroups)} RIBAP groups with the help of the pairwise ILPs.")
print(f"Of those new {len(assignedGroups)} groups, {coreGeneCounter} build core gene set including all {NUMSTRAINS} input genomes.")
