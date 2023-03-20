/*Comment section: */

process mmseqs2 {
  label 'mmseqs2'
  publishDir "${params.output}/03-mmseqs2", mode: 'copy', pattern: "mmseq2_result.csv" 
  publishDir "${params.output}/03-mmseqs2", mode: 'copy', pattern: "mmseq2_result_filtered.csv" 

  input: 
    file(fasta)

  output:
    path("mmseq2_result.csv")
    path("mmseq2_result_filtered.csv")

  script:
    """
    mkdir mmseq2
    cat *.faa > mmseq2/all_proteins.fa
    MMSEQDB="mmseq2/mmseq2.db"

    #Creating DB
    mmseqs createdb mmseq2/all_proteins.fa \$MMSEQDB
    #Creating Index#
    mmseqs createindex --threads ${task.cpus} \$MMSEQDB mmseq2/tmp
    #Starting MMSeqs2 Search
    mmseqs search --threads ${task.cpus} \$MMSEQDB \$MMSEQDB "\${MMSEQDB%.*}_result" mmseq2/tmp -a
    #Converting results
    mmseqs convertalis --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,qcov,tcov,evalue,bits" --threads ${task.cpus} \$MMSEQDB \$MMSEQDB "\${MMSEQDB%.*}_result" "\${MMSEQDB%.*}_result.csv"

    cat "\${MMSEQDB%.*}_result.csv" | awk '{if(\$3>0.6 && \$13>0.4){print \$0}}' > "\${MMSEQDB%.*}_result_filtered.csv"

    mv mmseq2/mmseq2_result.csv .
    mv mmseq2/mmseq2_result_filtered.csv .
    """
}

process mmseqs2tsv {
  label 'python3_high_mem'
  // publishDir "${params.output}/04-mmseqs2tsv", mode: 'copy', pattern: "*.tsv" 

  input: 
    file(mmseqs2)
    file(strain_ids)

  output:
    file("*.pkl")

  script:
    """
    #mkdir tsv
    mmseq2tsv.py ${mmseqs2} ${strain_ids} . ${params.chunks} #tsv 
    """
}

