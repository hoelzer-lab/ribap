import sys
import config
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import Bio.Alphabet
from Bio import SeqFeature as sf

"""This module contains two functions to rearrange a genome (fasta file) and, if given, a annotaion (Genbank file).
The function findPOI finds the breakpoint with respect to a certain gene product.
The function rearrange creates the rearranged genome and annoation file and saves is.

Author: Marie Lataretu
E-Mail: marie.lataretu@uni-jena.de
"""
def rearrange_fasta_file(genome_file, breakpoint, reverse_complement, output_name, contig_name=None):
    """This function rearranges and saves a genome (fasta file) and, if given, a annotation file (Genbank file)
    Required parameter: breakpoint (int), reverse_complement (True/False), output_name, genome_file, product_of_interest
    Optional parameter: contig name (in case of preliminary structural search in multiple contigs)
    """
    print ('start rearranging')

    ### read in the genome
    genome = None
    all_contigs_in_the_file = []
    with open(genome_file) as fas:
        for record in SeqIO.parse(fas, 'fasta'):
            if ('plasmid' not in record.description):
                if(contig_name):
                    if(contig_name == record.id):
                        genome = record
                else:
                    # no contig name given and assuming that there is only one contig
                    genome = record
                all_contigs_in_the_file.append(record)

    if(not genome):
        if(contig_name):
            sys.exit('end rearranging: no contig with the name ' + contig_name + ' found.')
    if(len(all_contigs_in_the_file) > 1 and not contig_name):
        sys.exit('end rearranging: genome consists of more than one non-plasmid contigs and it is not specified witch contig to rearrange')
    multiple_contigs = False
    if(len(all_contigs_in_the_file) > 1 and contig_name):
        multiple_contigs = True
        print('rearranging: multiple contigs')

    ### genome stuff
    if(reverse_complement):
        ## reverse complement and adjust breakpoint
        if(multiple_contigs):
            new_genome_sequence = genome.seq.reverse_complement()
            breakpoint = len(genome.seq) - breakpoint - 1  # annoying index things
            new_genome_sequence_1a = new_genome_sequence[breakpoint:]
            new_genome_sequence_1b = new_genome_sequence[:breakpoint]
        else:
            new_genome_sequence = genome.seq.reverse_complement()
            breakpoint = len(genome.seq) - breakpoint - 1   # annoying index things
            new_genome_sequence = new_genome_sequence[breakpoint:] + new_genome_sequence[:breakpoint]
    else:
        if(multiple_contigs):
            new_genome_sequence_1a = genome.seq[breakpoint:]
            new_genome_sequence_1b = genome.seq[:breakpoint]
        else:
            new_genome_sequence = genome.seq[breakpoint:] + genome.seq[:breakpoint]

    ## save the rearranges genome
    # if(contig_name):
    if(multiple_contigs):
        # new names
        genome.seq = new_genome_sequence_1a
        if(len(genome.id.split('_')) > 2):
            ## spades specific name
            print('spades specific name stuff')
            old_name_1a = genome.id.split('_')
            old_name_1a[1] = '1a'
            old_name_1a[3] = str(len(genome.seq))
            new_name_1a  = '_'.join(old_name_1a)

            old_name_1b = genome.id.split('_')
            old_name_1b[1] = '1b'
            old_name_1b[3] = str(len(new_genome_sequence_1b))
            new_name_1b = '_'.join(old_name_1b)

            genome.id = new_name_1a
            genome.description = new_name_1a
            genome.name = new_name_1a
            new_seq_record = SeqRecord(new_genome_sequence_1b, id=new_name_1b, name=new_name_1b,
                                       description=new_name_1b)

        elif(len(genome.id.split('.')) == 2):
            ## FLI specific name
            print('FLI specific name stuff')
            old_name_1a = genome.id.split('.')
            old_name_1a[0] = ''.join([old_name_1a[0], 'a'])
            new_name_1a = '.'.join(old_name_1a)

            new_des_1a = genome.description.split(' ')
            new_des_1a[0] = new_name_1a
            new_des_1a.append('length')
            new_des_1a.append(str(len(new_genome_sequence_1a)))

            old_name_1b = genome.id.split('.')
            old_name_1b[0] = ''.join([old_name_1b[0], 'b'])
            new_name_1b = '.'.join(old_name_1b)

            new_des_1b = genome.description.split(' ')
            new_des_1b[0] = new_name_1b
            new_des_1b.append('length')
            new_des_1b.append(str(len(new_genome_sequence_1b)))

            genome.id = new_name_1a
            genome.description = ' '.join(new_des_1a)
            genome.name = new_name_1a
            new_seq_record = SeqRecord(new_genome_sequence_1b, id=new_name_1b, name=new_name_1b,
                                       description=' '.join(new_des_1b))
        else:
            new_name_1a = '_'.join([genome.id, 'a', 'length', str(len(genome.seq))])
            new_name_1b = '_'.join([genome.id, 'b', 'length', str(len(new_genome_sequence_1b))])

            genome.id = new_name_1a
            genome.description = new_name_1a
            genome.name = new_name_1a
            new_seq_record = SeqRecord(new_genome_sequence_1b, id=new_name_1b, name=new_name_1b,
                                       description=new_name_1b)


        contig_index = None
        for i, record in enumerate(all_contigs_in_the_file):
            if(record.id == genome.id):
                contig_index = i

        all_contigs_in_the_file.insert(0, all_contigs_in_the_file.pop(contig_index))
        all_contigs_in_the_file.append(new_seq_record)
        SeqIO.write(all_contigs_in_the_file, config.OUTPUT_DIR + output_name + '_rearranged.fasta', 'fasta')
    else:
        genome.seq = new_genome_sequence

        SeqIO.write(genome, config.OUTPUT_DIR + output_name + '_rearranged.fasta', 'fasta')\

    print ('end rearranging')

def rearrange_genbank_file(anno_file, breakpoint, reverse_complement, output_name):
    print('start rearranging')
    if(anno_file):
        annotation = None
        ### read in the annotation file
        with open(anno_file) as gbk:
            for seq_record in SeqIO.parse(gbk, 'genbank'):
                if ('complete genome' in seq_record.description and ('plasmid' not in seq_record.description)):
                    annotation = seq_record

        if(not annotation):
            sys.exit('No complete genome in the Genbank file found.')

        for seq_feature in annotation.features:
            if (seq_feature.type == 'source'):
                ## get the genome length
                genome_length_anno = seq_feature.location.end
            else:
                ## discriminate the two different feature location classes (FeatureLocation and CompoundLocation)
                if(type(seq_feature.location) is sf.FeatureLocation):
                    if (breakpoint <= seq_feature.location.start):
                        # feature start behind breakpoint
                        new_start = seq_feature.location.start - breakpoint
                        new_end = seq_feature.location.end - breakpoint
                        seq_feature.location = sf.FeatureLocation(new_start, new_end, seq_feature.strand)
                    elif (breakpoint > seq_feature.location.end):
                        # feature end before breakpoint
                        new_start = genome_length_anno - breakpoint + seq_feature.location.start
                        new_end = genome_length_anno - breakpoint + seq_feature.location.end
                        seq_feature.location = sf.FeatureLocation(new_start, new_end, seq_feature.strand)
                    else:
                        # feature and breakpoint/product of interest are overlapping
                        new_start = seq_feature.end - breakpoint
                        new_end = genome_length_anno - (genome_length_anno - seq_feature.location.start)
                        seq_feature.location = sf.FeatureLocation(0, new_end, seq_feature.strand) + sf.FeatureLocation(
                        new_start, genome_length_anno, seq_feature.strand)
                elif(type(seq_feature.location) is sf.CompoundLocation):
                    # get the strand of the feature to calculate the new start correctly
                    if(seq_feature.location.strand is 1):
                        new_start = list(seq_feature.location)[0] - breakpoint
                    elif (seq_feature.location.strand is -1):
                        new_start = list(seq_feature.location)[-1] - breakpoint
                    else:
                        print('Mixed strands in one feature: ' + seq_feature)
                    new_end = new_start + len(seq_feature.location)
                    seq_feature.location = sf.FeatureLocation(new_start, new_end, seq_feature.strand)
                else:
                    sys.exit('end rearranging: something wired happened: a feature location is neither FeatureLocation nor CompoundLocation')

        if (reverse_complement):
            ## reverse complement and adjust breakpoint
            new_genome_sequence = annotation.seq.seq.reverse_complement()
            breakpoint = len(annotation.seq.seq) - breakpoint
            new_genome_sequence = new_genome_sequence[breakpoint:] + new_genome_sequence[:breakpoint]
        else:
            new_genome_sequence = annotation.seq[breakpoint:] + annotation.seq[:breakpoint]

        ## save the rearranges genome
        annotation.seq = new_genome_sequence
        # annotation.seq = Seq(str(annotation.seq), IUPAC.unambiguous_dna)

        ### write rearranged annotation
        SeqIO.write(annotation, config.OUTPUT_DIR + output_name + '_annotation_rear.gbk', 'genbank')

        print('end rearranging')

def rearrange_gff_file(gff_file, sequence_length, breakpoint, output_name):
    new_gff_file = '{}/{}'.format(config.OUTPUT_DIR, output_name)
    with open(gff_file, 'r') as f, open(new_gff_file, 'w') as f_out:
        reader = csv.reader(f, delimiter='\t')
        writer = csv.writer(f_out, delimiter='\t')
        for line in reader:
            if (len(line) is 0 or line[0].startswith('#')):
                writer.writerow(line)
            else:
                start = int(line[3])
                end = int(line[4])
                if (start >= breakpoint):
                    new_start = start - breakpoint + 1
                    new_end = end - breakpoint + 1
                    line[3] = new_start
                    line[4] = new_end
                    writer.writerow(line)
                elif (start < breakpoint and end >= breakpoint):
                    new_start_1 = sequence_length + start + 1 - breakpoint
                    new_end_1 = sequence_length
                    new_start_2 = 1
                    new_end_2 = 1 + end - breakpoint
                    line[0] = line[0] + '_1'
                    line[3] = new_start_1
                    line[4] = new_end_1
                    writer.writerow(line)
                    line[0] = line[0][:-1] + '2'
                    line[3] = new_start_2
                    line[4] = new_end_2
                    writer.writerow(line)
                else:
                    new_start = sequence_length - breakpoint + start
                    new_end = sequence_length - breakpoint + end
                    line[3] = new_start
                    line[4] = new_end
                    writer.writerow(line)

def get_sequence_length_fasta(single_fasta_file):
    with open(single_fasta_file) as f:
        record = SeqIO.read(f, 'fasta')
        return (len(record))

def get_sequence_length_genbank(genbank_file):
    with open(genbank_file) as gbk:
        for seq_record in SeqIO.parse(gbk, 'genbank'):
            if (('complete genome' in seq_record.description) and ('plasmid' not in seq_record.description)):
                annotation = seq_record

    if (not annotation):
        sys.exit('No complete genome in the Genbank file found.')
    for seq_feature in annotation.features:
        if (seq_feature.type == 'source'):
            ## get the genome length
            return (seq_feature.location.end)

if(__name__ == '__main__'):
    rearrange_fasta_file(17264, False, 'scaffold1_rearranged', '/mnt/dessertlocal/kons/Aspergillus_fumigatus_A1163/mitos/soapdenovo2_cap_sspace_reduced/0/sequence.fas-0', '')
