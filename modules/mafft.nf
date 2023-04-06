/*Comment section: */

process mafft {
  label 'mafft'
  publishDir "${params.output}/07-mafft", mode: 'copy', pattern: "*.aln" 

  input: 
    file(faa)

  output:
    file("*.aln")

  script:
    """
      for FAA in *.faa; do 
        mafft \${FAA} > "\$(basename \${FAA} .faa)"_mafft.aln
      done 
    """
}

