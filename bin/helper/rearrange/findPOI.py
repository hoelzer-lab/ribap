import sys
import os
import csv
import config
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature as sf

def find_poi_in_genbank_file(anno_with_POI_file, product_of_interest, save_sequence=False, output_name=None):
    non_plasmid_contig_list = []
    breakpoint = None
    reverse_complement = False

    ### search for the product of interest in the annotation file
    with open(anno_with_POI_file) as gbk:
        for seq_record in SeqIO.parse(gbk, 'genbank'):
            if ('complete genome' in seq_record.description):
                non_plasmid_contig_list.append(seq_record)
            elif('plasmid' not in seq_record.description):
                non_plasmid_contig_list.append(seq_record)

    if(len(non_plasmid_contig_list) != 0):
        for contig in non_plasmid_contig_list:
            contig_description = contig.description.replace(' ', '_')
            for seq_feature in contig.features:
                if ('product' in seq_feature.qualifiers):
                    if (product_of_interest in seq_feature.qualifiers['product'][0]):

                        breakpoint = seq_feature.location.start

                        ## check the strand
                        if (seq_feature.location.strand == -1):
                            reverse_complement = True

                        ##todo gleich reversekompl

                        if (save_sequence):
                            ## save the nucleotide sequence (e.g. to blast it against an another genome)
                            if (type(seq_feature.location) is sf.FeatureLocation):
                                ## file to save the nucleotide sequence of the product of interest
                                if (not os.path.exists(config.POI_OUTPUT_DIR)):
                                    os.makedirs(config.POI_OUTPUT_DIR)

                                if(output_name):
                                    POI_nc_file = config.POI_OUTPUT_DIR + output_name + '.fasta'
                                else:
                                    POI_nc_file = config.POI_OUTPUT_DIR + product_of_interest.replace(
                                        ' ', '_') + '_nc.fasta'

                                if(reverse_complement):
                                    new_description = '{}-{}-:{}-{}'.format(contig_description, 'revcomp',
                                                                            str(seq_feature.location.start + 1),
                                                                            str(seq_feature.location.end))

                                    new_nc_record = SeqRecord(
                                        contig.seq[seq_feature.location.start:seq_feature.location.end].reverse_complement(),
                                        id=product_of_interest, description=new_description)
                                else:
                                    new_description = '{}:{}-{}'.format(contig_description,
                                                                        str(seq_feature.location.start + 1),
                                                                        str(seq_feature.location.end))
                                    new_nc_record = SeqRecord(
                                        contig.seq[seq_feature.location.start:seq_feature.location.end],
                                         id=product_of_interest, description=new_description)

                                SeqIO.write(new_nc_record, POI_nc_file, 'fasta')
                            else:
                                print ('product of interest has a CompoundLocation; this is not implemented yet')
                                sys.exit(0)
    else:
        sys.exit('No complete genome in the Genbank file found.')
    if(breakpoint is None):
        sys.exit('Product of interesst not found.')
    return ([int(breakpoint), reverse_complement])

def find_poi_in_fasta_file(genome_fasta_file, seqeunce_of_interest_fasta_file, output_blast_name):
    reverse_complement = False
    save_contig = False

    ### read in the genome file
    num_of_records = 0
    with open(genome_fasta_file) as fas:
        for record in SeqIO.parse(fas, 'fasta'):
            if ('plasmid' not in record.description):
                genome = record
                num_of_records = num_of_records + 1
    if (num_of_records > 1):
        save_contig = True

    if (not os.path.exists(config.BLAST_OUTPUT_DIR)):
        os.makedirs(config.BLAST_OUTPUT_DIR)
    # out_blast_file = config.BLAST_OUTPUT_DIR + strain_name_with_POI + '_poi_vs_' + strain_name_to_search_in + '_genome.tab'
    out_blast_file = config.BLAST_OUTPUT_DIR + output_blast_name +'.tab'

    if (not (os.path.isfile(genome_fasta_file + '.nhr') and os.path.isfile(genome_fasta_file + '.nsq') and os.path.isfile(genome_fasta_file + '.nin'))):
        print ('start makeblastdb')
        makeblastdbStr = config.BLAST_BIN_PATH + 'makeblastdb -dbtype nucl -in ' + genome_fasta_file
        os.system(makeblastdbStr)
        print ('end makeblastdb')

    print ('start blast')
    blastStr = config.BLAST_BIN_PATH + 'blastn -out ' + out_blast_file + ' -outfmt \'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\' -query ' + seqeunce_of_interest_fasta_file + ' -db ' + genome_fasta_file + ' -evalue 1e-10'
    os.system(blastStr)
    print ('end blast')

    ## read in the blast result
    with open(out_blast_file) as csvf:
        table = [row for row in csv.reader(csvf, delimiter='\t')]
    if (len(table) is 1):

        ## check the stand of the product of interest
        if (int(table[0][8]) > int(table[0][9])):
            reverse_complement = True
        # return the breakpoint, if is is reverse complemented and on which contig the breakpoint is, if necessery
        if(save_contig):
            return ([int(table[0][8])-1, reverse_complement, table[0][1]])
        else:
            return ([int(table[0][8])-1, reverse_complement])
    elif(len(table) is 0):
        print ('no blast hit')
        sys.exit(0)
    else:
        print ('end find product of interest: found more than one blast hit; not implemented yet')
        sys.exit(0)

# def findPOI(strain_name_with_POI, genome_with_POI_file, anno_with_POI_file, product_of_interest, strain_name_to_search_in=None, genome_to_search_in=None):
#     """This function returns the breakpoint with respect to a certain gene product of interest (POI)
#     and if the nucleotide sequence must be reverse complemented.
#     Using the required parameters the function searches the annotation file by name. Using additionally the
#     optional parameters the function performs a structural search in the genome_to_search_in using blast and
#     genome_with_POI_file and anno_with_POI_file as reference.
#     Required parameters: strain_name_with_POI, genome_with_POI_file, anno_with_POI_file, product_of_interest
#     Optional parameters: strain_name_to_search_in, genome_to_search_in
#     """
#     print ('start find product of interest')
#
#     breakpoint = None
#     reverse_complement = False
#
#     ### read in the genome with the product of interest
#     num_of_records = 0
#     with open(genome_with_POI_file) as fas:
#         for record in SeqIO.parse(fas, 'fasta'):
#             if ('plasmid' not in record.description):
#                 genome_with_POI = record
#                 num_of_records = num_of_records  + 1
#     if(num_of_records > 1):
#         sys.exit('end find product of interest: genome consists of more than one non-plasmide contig; not implemented yet')
#
#     ### search for the product of interest in the annotation file
#     with open(anno_with_POI_file) as gbk:
#         for seq_record in SeqIO.parse(gbk, 'genbank'):
#             if ('complete genome' in seq_record.description):
#                 for seq_feature in seq_record.features:
#                     if ('product' in seq_feature.qualifiers):
#                         if (product_of_interest in seq_feature.qualifiers['product'][0]):
#
#                             breakpoint = seq_feature.location.start
#
#                             ## check the strand
#                             if(seq_feature.location.strand == -1):
#                                 reverse_complement = True
#
#                             if(genome_to_search_in):
#                                 ## save the nucleotide sequence to blast it against another genome
#                                 if(type(seq_feature.location) is sf.FeatureLocation):
#                                     ## file to save the nucleotide sequence of the product of interest
#                                     if (not os.path.exists(config.POI_OUTPUT_DIR)):
#                                         os.makedirs(config.POI_OUTPUT_DIR)
#                                     POI_nc_file = config.POI_OUTPUT_DIR + strain_name_with_POI + '_' + product_of_interest.replace(
#                                         ' ', '_') + '_nc.fasta'
#                                     new_nc_record = SeqRecord(genome_with_POI.seq[seq_feature.location.start:seq_feature.location.end], id=product_of_interest)
#                                     SeqIO.write(new_nc_record, POI_nc_file, 'fasta')
#                                 else:
#                                     print ('product of interest has a CompoundLocation; this is not implemented yet')
#                                     sys.exit(0)
#
#     if(breakpoint and genome_to_search_in):
#         ### product of interest found in the annotation; now blast the nucleotide sequence against the second genome
#
#         if (not os.path.exists(config.BLAST_OUTPUT_DIR)):
#             os.makedirs(config.BLAST_OUTPUT_DIR)
#         out_blast_file = config.BLAST_OUTPUT_DIR + strain_name_with_POI + '_poi_vs_' + strain_name_to_search_in + '_genome.tab'
#
#         ## prepare blast database if necessary
#         if (not (os.path.isfile(genome_to_search_in + '.nhr') and os.path.isfile(genome_to_search_in + '.nsq') and os.path.isfile(genome_to_search_in + '.nin'))):
#             print ('start makeblastdb')
#             makeblastdbStr = config.BLAST_BIN_PATH + 'makeblastdb -dbtype nucl -in ' + genome_to_search_in
#             os.system(makeblastdbStr)
#             print ('end makeblastdb')
#
#         ## blast product of interest found against the second genome
#         print ('start blast')
#         blastStr = config.BLAST_BIN_PATH + 'blastn -out ' + out_blast_file + ' -outfmt \'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\' -query ' + POI_nc_file + ' -db ' + genome_to_search_in + ' -evalue 1e-10'
#         os.system(blastStr)
#         print ('end blast')
#
#         ## read in the blast result
#         with open (out_blast_file) as csvf:
#             table = [row for row in csv.reader(csvf, delimiter='\t')]
#
#         if(len(table) is 1):
#             print ('end find product of interest')
#
#             ## check the stand of the product of interest
#             if(int(table[0][8]) > int(table[0][9])):
#                 reverse_complement = True
#             # return the breakpoint, if is is reverse complemented and on which contig the breakpoint is
#             return([int(table[0][8]), reverse_complement, table[0][1]])
#         else:
#             print ('end find product of interest: found more than one blast hit; not implemented yet')
#             sys.exit(0)
#
#     elif(breakpoint and genome_to_search_in is None):
#         ### one genome; product of interest found
#         print ('end find product of interest')
#         return ([int(breakpoint), reverse_complement])
#
#     else:
#         ### product of interest not found in the annotation
#         print ('end find product of interest: product of interest not found')
#         sys.exit(0)

