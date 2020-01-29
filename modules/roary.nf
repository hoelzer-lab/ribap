/*Comment section: */

process roary {
  label 'roary'
  publishDir "${params.output}/roary", mode: 'copy', pattern: "${ident}" 

  input: 
    tuple val(ident), file(gff)

  output:
    tuple val(ident), file(ident)

  script:
    """
    roary -e --mafft -p ${task.cpus} -v -i ${ident} -r *.gff -f "${ident}" &> ribap_roary_"${ident}".log
    """
}

