/*Comment section: */

process prepare_msa {
  label 'python3'
  errorStrategy{task.exitStatus=1 ?'ignore':'terminate'}
  publishDir "${params.output}/07-msa/", mode: 'copy', pattern: "msa/*.faa" 

  input: 
    tuple val(ident), path(holy_ribap_csv)
    path(faa)

  output:
    file("msa/*.faa")

  script:
    """
      mkdir faa
      cp *.faa faa/

      mkdir msa
      create_msa_tree.py . ${holy_ribap_csv}
      #mv msa/*.faa .
    """
}

