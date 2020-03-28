#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
nice docstring for this script
"""

from __future__ import print_function
__author__ = "Kevin Lamkiewicz"

import sys
import os
from collections import defaultdict
import colorsys
import random
import glob


projectDir = sys.argv[1]
ribapFile = f'{projectDir}/holy_python_ribap_95.csv'
ribapHash = defaultdict(list)

column2strain = {}

with open(ribapFile, 'r') as inputStream:
    while 1:
        currentLine = inputStream.readline().rstrip("\n")
        if not currentLine:
            break
        
        currentArray = currentLine.split("\t")
        # do this once at the header line
        if currentLine.startswith("Cluster"):
            for idx, strain in enumerate(currentArray[3:]):
                column2strain[idx] = strain
            continue

        # do this block for all line EXCEPT the header line
        ribapGroup = currentArray[0]
        geneDesc = currentArray[1]
        geneName = currentArray[2]
        genes = currentArray[3:]
        geneTable = {}

        for idx, geneID in enumerate(genes):
            strain = column2strain[idx]
            geneTable[strain] = geneID
        ribapHash[ribapGroup] = [geneName, geneDesc, geneTable]

individualAnnotation = f'{projectDir}/ribap_individual_annotation.csv'
individualHash = defaultdict(dict)

with open(individualAnnotation, 'r') as inputStream:
    while 1:
        currentLine = inputStream.readline().rstrip("\n")
        if not currentLine:
            break
        
        if currentLine.startswith("Cluster"):
            continue

        currentArray = currentLine.split("\t")
        ribapGroup = currentArray[0]
        strains = currentArray[1:]

        individualHash[ribapGroup] = {column2strain[idx]: annotation for idx,annotation in enumerate(strains)}

strains = list(column2strain.values())
NUMSTRAINS = len(strains)
POS = list(range(14,14+NUMSTRAINS))

IDENT = []
roaryCluster = {}
for roaryFile in glob.glob(f'{projectDir}/roary/*/gene_presence_absence.csv'):
    paralogs = []
    roarySim = os.path.basename(os.path.dirname(roaryFile))
    roaryCluster[roarySim] = {}
    IDENT.append(roarySim)
    currentRoaryTable = roaryCluster[roarySim]
    with open(roaryFile, 'r') as inputStream:
        for clusterID, line in enumerate(inputStream):
            if line.startswith('"Gene'):
                continue
            
            array = line.rstrip('\n').split('","')
            roarygroup = []
            for column in POS:
                if array[column] and not array[column] == '"':
                    if '\t' in array[column]:
                        array[column], paralogs = array[column].split('\t')[0], array[column].split('\t')[1:]
                        currentRoaryTable[array[column].rstrip('"')] = clusterID
                        for idx, para in enumerate(paralogs):
                            currentRoaryTable[para.rstrip('"')] = f'{clusterID}.{idx}'
                    else:
                        currentRoaryTable[array[column].rstrip('"')] = clusterID

randomColors = set()
while len(randomColors) != 5*NUMSTRAINS:
    randomColors.add(colorsys.hsv_to_rgb(random.random(), random.random(), random.random()))

HTMLHEADER = """
<!doctype html>
<html lang="en">

<head>
    <meta http-equiv="Content-type" content="text/html; charset=utf-8">
    <meta name="viewport" content="width=device-width,initial-scale=1,user-scalable=no">
    <title>RIBAP</title>
    <link rel="stylesheet" type="text/css" href="jquery.dataTables.min.css">

    <style type="text/css" media="screen">
        @import url('//cdn.datatables.net/1.10.2/css/jquery.dataTables.css');
        td.details-control {
            background: url('http://www.datatables.net/examples/resources/details_open.png') no-repeat center center;
            cursor: pointer;
        }

        tr.shown td.details-control {
            background: url('http://www.datatables.net/examples/resources/details_close.png') no-repeat center center;
        }
    </style>

    <script type="text/javascript" language="javascript" src="jquery-3.3.1.min.js"></script>
    <script type="text/javascript" language="javascript" src="jquery.dataTables.min.js"></script>
    <script type="text/javascript" class="init">

        function format(value) {
            return '<tr>' + value + '</tr>';
        }
        $(document).ready(function () {
            var table = $('#ribap').DataTable({});

            // Add event listener for opening and closing details
            $('#ribap').on('click', 'td.details-control', function () {
                var tr = $(this).closest('tr');
                var row = table.row(tr);

                if (row.child.isShown()) {
                    // This row is already open - close it
                    row.child.hide();
                    tr.removeClass('shown');
                } else {
                    // Open this row
                    row.child(format(tr.data('child-value'))).show();
                    tr.addClass('shown');
                }
            });
        });

    </script>
</head>
"""
            
HTMLBODY = """
<body>

    <table id="ribap" class="display nowrap" cellspacing="0" width="100%">
        <thead>
            <tr>
                <th></th>
                <th>Group ID</th>    
                <th>Gene Name</th>
                <th>Gene Description</th>
                <th>Group Size</th>
"""



for strain in strains:
    HTMLBODY += f"<th>{strain}</th>\n"

HTMLBODY += """
    </tr>
</thead>
<tbody>
"""

for group, content in ribapHash.items():
    

    uniqueStrains = set(content[2].values())
    if "NA" in uniqueStrains:
        uniqueStrains.remove("NA")
    
    if not uniqueStrains:
        continue
        print(group)
        exit(0)
    
    roary60 = { gene : group for gene, group in roaryCluster['60'].items() if gene in uniqueStrains }
    roary70 = { gene : group for gene, group in roaryCluster['70'].items() if gene in uniqueStrains }
    roary80 = { gene : group for gene, group in roaryCluster['80'].items() if gene in uniqueStrains }
    roary90 = { gene : group for gene, group in roaryCluster['90'].items() if gene in uniqueStrains }
    roary95 = { gene : group for gene, group in roaryCluster['95'].items() if gene in uniqueStrains }
    
    combinedSet = list(roary60.values()) + list(roary70.values()) + list(roary80.values()) + list(roary90.values()) + list(roary95.values())
    sampledColors = random.sample(randomColors, len((combinedSet)))
    
    colors60 = dict(zip(set(roary60.values()), sampledColors))
    for usedColor in colors60.values():
        sampledColors.remove(usedColor)
    sampledColors60 = {gene: colors60[group] for gene, group in roary60.items()}
    sampledColors60.update({'NA' : (1,1,1)})

    colors70 = dict(zip(set(roary70.values()), sampledColors))
    for usedColor in colors70.values():
        sampledColors.remove(usedColor)
    sampledColors70 = {gene: colors70[group] for gene, group in roary70.items()}
    sampledColors70.update({'NA' : (1,1,1)})

    colors80 = dict(zip(set(roary80.values()), sampledColors))
    for usedColor in colors80.values():
        sampledColors.remove(usedColor)
    sampledColors80 = {gene: colors80[group] for gene, group in roary80.items()}
    sampledColors80.update({'NA' : (1,1,1)})

    colors90 = dict(zip(set(roary90.values()), sampledColors))
    for usedColor in colors90.values():
        sampledColors.remove(usedColor)
    sampledColors90 = {gene: colors90[group] for gene, group in roary90.items()}
    sampledColors90.update({'NA' : (1,1,1)})

    colors95 = dict(zip(set(roary95.values()), sampledColors))
    for usedColor in colors95.values():
        sampledColors.remove(usedColor)
    sampledColors95 = {gene: colors95[group] for gene, group in roary95.items()}
    sampledColors95.update({'NA' : (1,1,1)})
    
    
        
    if len(uniqueStrains) <= 2:
        tree = ""
    else:
        #tree =f"<img src={projectDir}/tree/{group}_mafft_tree.svg>"
        tree =f"<img src=../tree/{group}_mafft_tree.svg>"


    HTMLBODY += '<tr data-child-value="<th>Strain</th><th>Gene Name</th><th>Gene Description</th><th>Roary 60</th><th>Roary 70</th><th>Roary 80</th><th>Roary 90</th><th>Roary 95</th><th>Phylogenetic Tree</th>'
    i = 0
    for strain, individualAnno in individualHash[group].items():
        if individualAnno.split(" // ")[0] == '--' and individualAnno.split(" // ")[1] == "NA":
            continue
        if i == 0:
            HTMLBODY += f'<tr><td>{strain}</td><td>{individualAnno.split(" // ")[0]}</td><td>{individualAnno.split(" // ")[1]}</td>\n<td bgcolor=\'rgb{sampledColors60[content[2][strain]]}\'> </td><td bgcolor=\'rgb{sampledColors70[content[2][strain]]}\'> </td><td bgcolor=\'rgb{sampledColors80[content[2][strain]]}\'> </td><td bgcolor=\'rgb{sampledColors90[content[2][strain]]}\'> </td><td bgcolor=\'rgb{sampledColors95[content[2][strain]]}\'> </td>\n<td rowspan=\'{len(individualHash[group].keys())}\'>{tree}<p align=\'center\'><a href=\'../msa/{group}_mafft.aln\'>Multiple Sequence Alignment</a><br><a href=\'../msa/{group}_mafft_tree.nwk\'>Newick Tree</a></p></td></tr>\n'
            i += 1
        else:
            HTMLBODY += f'<tr><td>{strain}</td><td>{individualAnno.split(" // ")[0]}</td><td>{individualAnno.split(" // ")[1]}</td><td bgcolor=\'rgb{sampledColors60[content[2][strain]]}\'> </td><td bgcolor=\'rgb{sampledColors70[content[2][strain]]}\'> </td><td bgcolor=\'rgb{sampledColors80[content[2][strain]]}\'> </td><td bgcolor=\'rgb{sampledColors90[content[2][strain]]}\'> </td><td bgcolor=\'rgb{sampledColors95[content[2][strain]]}\'> </td></tr>\n'
    
    HTMLBODY += '">\n<td class="details-control"></td>\n'



    HTMLBODY += f'<td>{group}</td><td>{content[0]}</td><td>{content[1]}</td><td>{len(uniqueStrains)}</td>'
    
    for strain, geneID in content[2].items():
        HTMLBODY += f'<td>{geneID}</td>\n'
    
    HTMLBODY += '</tr>'

HTMLBODY += """
        </tbody>
    </table>

</body>
</html>
"""

print(HTMLHEADER)
print(HTMLBODY)