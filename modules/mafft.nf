/*Comment section: */

process mafft {
  label 'mafft'
  publishDir "${params.output}/mafft", mode: 'copy', pattern: "*.aln" 

  input: 
    file(faa)

  output:
    file("*.aln")

  script:
    """
      mafft ${faa} > "\$(basename ${faa} .faa)"_mafft.aln 
    """
}

