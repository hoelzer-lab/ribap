/*Comment section: */

process fasttree {
  label 'fasttree'
  publishDir "${params.output}/fasttree", mode: 'copy', pattern: "*.nwk" 

  input: 
    file(aln)

  output:
    file("*.nwk")

  script:
    """
    fasttree *.aln > "\$(basename ${aln} _mafft.aln)"_tree.nwk 
    """
}

