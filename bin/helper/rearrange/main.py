import sys
import argparse
from findPOI import find_poi_in_genbank_file
from findPOI import find_poi_in_fasta_file
from rearranging import rearrange_fasta_file
from rearranging import rearrange_genbank_file
from rearranging import rearrange_gff_file
from rearranging import get_sequence_length_fasta
from rearranging import get_sequence_length_genbank

"""This script rearranges a genome (fasta file) and, if given, a annotaion (Genbank file) with respect to a product
of interest (POI).

Output paths can be specified in config.py. 

Usage:
python main.py "strain_name_with_POI" "genome_with_POI_file" "anno_with_POI_file" "product_of_interest"
["strain_name_to_search_in"] ["genome_to_search_in"] ["anno_to_search_in"]

Examples:
a) given a genome and annotation file of a species called 'fuubarus' and a annotated product 'muh' 
-> python main.py "fuubarus" "fuubarus.fasta" "fuubarus.gbf" "muh"

b) given a genome file of a species called 'barbarus' and a product 'muh', that is not annotated. take the genome and
annotation of 'fuubarus' to get the sequence of 'muh' and run blast to find it in 'barbarus'
-> python main.py "fuubarus" "fuubarus.fasta" "fuubarus.gbf" "muh" "barbarus" "barbarus.fasta"

c) given a genome file of a species called 'barbarus' and a product 'muh', that is not annotated namely in the 
annotation. 
take the genome and annotation of 'fuubarus' to get the sequence of 'muh' and run blast to find it in 'barbarus'.
-> python main.py "fuubarus" "fuubarus.fasta" "fuubarus.gbf" "muh" "barbarus" "barbarus.fasta" "barbarus.gbf"

Used tools: ncbi-blast-2.4.0+


Author: Marie Lataretu
E-Mail: marie.lataretu@uni-jena.de
"""

parser = argparse.ArgumentParser(description='command line tool to rearrange fasta and genbank files. and in some future version also gff files')

parser.add_argument('-structuralSearch', action='store_true', help='do a structural search to find the breakpoint')
parser.add_argument('-nameSearch', action='store_true', help='find the breakpoint via the name of a product in a genbank file')
parser.add_argument('-rearrangeFasta', action='store_true', help='rearrange a fasta file')
parser.add_argument('-rearrangeGenbank', action='store_true', help='rearrange a genbank file')
parser.add_argument('-rearrangeGff', action='store_true', help='rearrange a gff file')


parser.add_argument('-outputName', type=str, metavar='STR', help='name of the output file')
parser.add_argument('-genomeFile', type=str, metavar='FILE', help='fasta file of the genome')
parser.add_argument('-annotationFile', type=str, metavar='FILE', help='genbank file of the genome')
parser.add_argument('-gffFile', type=str, metavar='FILE', help='gff file of the genome')

parser.add_argument('-POIName', type=str, metavar='STR', help='name of the product of interest')
parser.add_argument('-sequenceOfInterest', type=str, metavar='FILE', help='nucleotide sequence for structural search')
parser.add_argument('-blastOutName', type=str, metavar='STR', help='name for the outputfile of the structural search')

parser.add_argument('-POIOutName', type=str, metavar='STR', help='name for the sequence outputfile of product of interest')
parser.add_argument('-savePOISequence', action='store_true', help='save the sequence of the name search')

# args = parser.parse_args(['-nameSearch', '-annotationFile', 'TEST_Cavium_annotation_small.gbk', '-POIName', '\'delta-aminolevulinic acid dehydratase\''])
# args = parser.parse_args(['-nameSearch', '-annotationFile', '/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_annotation.gbk',
#                           '-POIName', '\'delta-aminolevulinic acid dehydratase\'', '-savePOISequence', '-POIOutName', 'testseqpsittaci'])
# args = parser.parse_args(['-structuralSearch', '-genomeFile', '/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_genome.fasta',
#                           '-sequenceOfInterest', '/mnt/prostlocal/marie/chlamydia_comparison/test/POI/testseqpsittaci.fasta', '-blastOutName', 'blasttestpsittaci'])
# args = parser.parse_args(['-nameSearch', '-annotationFile', '/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_annotation.gbk',
#                           '-POIName', '\'delta-aminolevulinic acid dehydratase\'',
#                             '-rearrangeFasta', '-genomeFile', '/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_genome.fasta',
#                           '-outputName', 'testgenome_name_psittaci'])
# args = parser.parse_args(['-structuralSearch', '-sequenceOfInterest', '/mnt/prostlocal2/projects/chlamydia_comparative_study/lisa/rearrange/POI/Cga_hemB.fasta',
#                           '-blastOutName', 'blasttest_gall',
#                             '-rearrangeFasta', '-genomeFile', '/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_gallinacea/14DC101_Cgallinacea_nodes1+2+3.fasta',
#                           '-outputName', 'testgenome_struc_gall_multiple_stuff'])
# args = parser.parse_args(['-nameSearch', '-annotationFile', '/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_annotation.gbk',
#                           '-POIName', '\'delta-aminolevulinic acid dehydratase\'',
#                           '-rearrangeGenbank', '-outputName', 'testgenbank_name_psittaci'])
# args = parser.parse_args(['-nameSearch', '-annotationFile', '/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_annotation.gbk',
#                           '-POIName', '\'delta-aminolevulinic acid dehydratase\'',
#                           '-rearrangeGff', '-gffFile', '/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_annotation.gff', '-outputName', 'testgff_name_psittaci'])
# args = parser.parse_args([])
# args = parser.parse_args(['-h'])
args = parser.parse_args()
# print (args)

searchResult = None

if(args.structuralSearch):
    if(args.nameSearch):
        sys.exit('You can perform either structural search or a name search.')
    if(not(args.genomeFile or args.sequenceOfInterest or args.blastOutName)):
        sys.exit('For a structural search you need -genomeFile, -sequenceOfInterest and -blastOutName.')

    searchResult = find_poi_in_fasta_file(args.genomeFile, args.sequenceOfInterest, args.blastOutName)

if(args.nameSearch):
    # optional POIOutName, savePOISequence
    if(args.structuralSearch):
        sys.exit('You can perform either structural search or a name search.')
    if(not(args.annotationFile or args.POIName)):
        sys.exit('For a name search you need -annotationFile and -product_of_interest.')
    if(args.POIOutName and not args.savePOISequence):
        sys.exit('-POIOutName makes only sense with -savePOISequence')

    searchResult = find_poi_in_genbank_file(args.annotationFile, args.POIName.replace('\'', ''), args.savePOISequence, args.POIOutName)

if(args.rearrangeFasta):
    if(not args.structuralSearch):
        sys.exit('You need to search first with -structuralSearch.')
    if(not(args.genomeFile or args.outputName)):
        sys.exit('You need to define a fasta file with -genomeFile and a name for the output with -outputName.')
    if(searchResult):
        if(len(searchResult) is 3):
            rearrange_fasta_file(args.genomeFile, searchResult[0], searchResult[1], args.outputName, searchResult[2])
        else:
            rearrange_fasta_file(args.genomeFile, searchResult[0], searchResult[1], args.outputName)
    else:
        sys.exit('Something went wrong with the search.')

if(args.rearrangeGenbank):
    if(not args.nameSearch):
        sys.exit('You need to search first with -nameSearch.')
    if(not(args.annotationFile or args.outputName)):
        sys.exit('You need to define a genbank file with -annotationFile and a name for the output with -outputName.')

    if(searchResult):
        rearrange_genbank_file(args.annotationFile, searchResult[0], searchResult[1], args.outputName)
    else:
        sys.exit('Something went wrong with the search.')

if(args.rearrangeGff):
    if(not (args.nameSearch or args.structuralSearch)):
        sys.exit('You need to specify a type of search with either -nameSearch or -structuralSearch.')
    if(not (args.gffFile or args.outputName)):
        sys.exit('You need to define a gff file with -gffFile and a name for the output with -outputName.')

    if (searchResult):
        if(args.genomeFile):
            rearrange_gff_file(args.gffFile, get_sequence_length_fasta(args.genomeFile), searchResult[0], args.outputName)
        else:
            rearrange_gff_file(args.gffFile, get_sequence_length_genbank(args.annotationFile), searchResult[0], args.outputName)
    else:
        sys.exit('Something went wrong with the search.')

# def main(argv):
#     ### check and pares input parameters
#     if(len(argv) == 4 or len(argv) == 6 or len(argv) == 7):
#         strain_name_with_POI = argv[0].replace(' ', '_')
#         genome_with_POI_file = argv[1]
#         anno_with_POI_file = argv[2]
#         product_of_interest = argv[3]
#
#     else:
#         sys.exit('4, 6 or 7 arguments are expected')
#
#
#     if(len(argv) == 4):
#         ### search by name
#         breakpoint_and_reverse_complement = findPOI(strain_name_with_POI, genome_with_POI_file, anno_with_POI_file, product_of_interest)
#         rearrange_fasta_file(breakpoint_and_reverse_complement[0], breakpoint_and_reverse_complement[1], strain_name_with_POI, genome_with_POI_file, product_of_interest, anno_with_POI_file)
#
#     elif(len(argv) > 4):
#         ### structural search
#         strain_name_to_search_in = argv[4].replace(' ', '_')
#         genome_to_search_in = argv[5]
#
#         breakpoint_reverseComplement_contig = findPOI(strain_name_with_POI, genome_with_POI_file, anno_with_POI_file, product_of_interest,
#                              strain_name_to_search_in, genome_to_search_in)
#
#         if(len(argv) == 6):
#             ## no second annotation file given
#             rearrange_fasta_file(breakpoint_reverseComplement_contig[0], breakpoint_reverseComplement_contig[1], strain_name_to_search_in,
#                                  genome_to_search_in, product_of_interest, contig_name=breakpoint_reverseComplement_contig[2])
#
#         else:
#             ## rearrange the second annotation file
#             anno_to_search_in = argv[6]
#             rearrange_fasta_file(breakpoint_reverseComplement_contig[0], breakpoint_reverseComplement_contig[1], strain_name_to_search_in,
#                                  genome_to_search_in, product_of_interest, breakpoint_reverseComplement_contig[2])
#             rearrang_genbank_file(breakpoint_reverseComplement_contig[0], breakpoint_reverseComplement_contig[1], anno_to_search_in, strain_name_to_search_in)

# if(__name__ == '__main__'):
    # print (sys.argv)
    # main(sys.argv[1:])
    # para = ["6BC_Cpsittaci",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_genome.fasta",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_annotation_edit.gbk",
    #         "delta-aminolevulinic acid dehydratase"]

    # para = ["10DC88_Cavium",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_avium/10DC88_Cavium_genome.fasta",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_avium/10DC88_Cavium_annotation.gbk",
    #         "delta-aminolevulinic acid dehydratase",
    #         "11DC096_Cavium",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_avium/11DC096_Cavium_chromosome_node1.fasta"]

    # para = ["6BC_Cpsittaci",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_genome.fasta",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_annotation_edit.gbk",
    #         "delta-aminolevulinic acid dehydratase",
    #         "02DC15_Cpsittaci",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/02DC15_Cpsittaci_genome.fasta",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/02DC15_Cpsittaci_annotation.gbk"]
    #
    # main(para)
