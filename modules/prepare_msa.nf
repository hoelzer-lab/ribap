/*Comment section: */

process prepare_msa {
  label 'python3'
  errorStrategy{task.exitStatus=1 ?'ignore':'terminate'}
  publishDir "${params.output}/msa/", mode: 'copy', pattern: "*.faa" 

  input: 
    tuple val(ident), file(holy_ribap_csv)
    file(faa)

  output:
    file("*.faa")

  script:
    """
      mkdir msa
      create_msa_tree.py . ${holy_ribap_csv}
      mv msa/*.faa .
    """
}

