/*Comment section: */

process fasttree {
  label 'fasttree'
  publishDir "${params.output}/09-fasttree", mode: 'copy', pattern: "*.nwk" 

  input: 
    file(aln)

  output:
    file("*.nwk")

  script:
    """
    for ALN in *.aln; do
      fasttree \${ALN} > "\$(basename \${ALN} _mafft.aln)"_tree.nwk 
    done
    """
}

