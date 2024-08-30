/*Comment section: */

process rename {
  label 'basics'
  publishDir "${params.output}/00-rename", mode: 'copy', pattern: "*_RENAMED.fasta" 
  tag "$name"

  input: 
    tuple val(name), file(fasta)

  output:
    tuple val(name), val("${name}_RENAMED"), file("*_RENAMED.fasta")

  script:
    """
    rename.sh ${fasta}
    """
}

