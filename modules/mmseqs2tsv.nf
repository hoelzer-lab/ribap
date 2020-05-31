/*Comment section: */

process mmseq2tsv {
  label 'python3'
  publishDir "${params.output}/blast2tsv", mode: 'copy', pattern: "*.tsv" 

  input: 
    file(mmseqs2)
    file(strain_ids)

  output:
    file("*.tsv")

  script:
    """
    #mkdir tsv
    blast2tsv.py ${mmseqs2} ${strain_ids} . #tsv 
    """
}

