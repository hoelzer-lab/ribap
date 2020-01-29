/*Comment section: */

process mmseqs2 {
  label 'mmseqs2'
  publishDir "${params.output}/mmseqs2", mode: 'copy', pattern: "mmseq2_result.csv" 

  input: 
    file(fasta)

  output:
    file("mmseq2_result.csv")

  script:
    """
    mkdir mmseq2
    cat *.faa > mmseq2/all_proteins.fa
    MMSEQDB="mmseq2/mmseq2.db"

    #Creating DB
    mmseqs createdb mmseq2/all_proteins.fa \$MMSEQDB
    #Creating Index#
    mmseqs createindex \$MMSEQDB mmseq2/tmp
    #Starting MMSeqs2 Search
    mmseqs search \$MMSEQDB \$MMSEQDB "\${MMSEQDB%.*}_result" mmseq2/tmp -a
    #Converting results
    mmseqs convertalis "\$MMSEQDB" "\$MMSEQDB" "\${MMSEQDB%.*}_result" "\${MMSEQDB%.*}_result.csv"

    #create single .tsv
    MMSEQ="\${MMSEQDB%.*}_result.csv"
    mv mmseq2/mmseq2_result.csv .
    """
}

