/*Comment section: */

process filter_alignment {
  label 'filter_alignment'
  publishDir "${params.output}/msa", mode: 'copy', pattern: "*core.aln" 

  input: 
    file(aln)
    file(strain_ids)

  output:
    file("*core.aln") optional true

  script:
    """
    NUM=\$(grep -c '>' ${aln})
    STRAINS=\$(cat ${strain_ids} | wc -l)
    if [ \$NUM -eq \$STRAINS ]; then
      cp ${aln} "\$(basename ${aln} .aln)"_core.aln
    fi
    """
}

