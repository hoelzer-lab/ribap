/*Comment section: */

process rename {
  label 'basics'
  publishDir "${params.output}/rename", mode: 'copy', pattern: "${name}_RENAMED.fasta" 

  input: 
    tuple val(name), file(fasta)

  output:
    tuple val("${name}_RENAMED"), file("${name}_RENAMED.fasta")

  script:
    """
    rename.sh ${fasta}
    """
}

