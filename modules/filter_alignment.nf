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
    for file in ${aln}; do
      NUM=\$(grep -c '>' \$file)
      STRAINS=\$(cat ${strain_ids} | wc -l)
      if [ \$NUM -eq \$STRAINS ]; then
        cp \${file} "\$(basename \${file} .aln)"_core.aln
        sed -r -i '/>/ s/_[^_]+\$//' "\$(basename \${file} .aln)"_core.aln
      fi
    done
    """
}

