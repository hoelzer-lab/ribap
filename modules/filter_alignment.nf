/*Comment section: */

process filter_alignment {
  label 'filter_alignment'
  publishDir "${params.output}/08-msa", pattern: "FINAL*core.aln" 

  input: 
    file(aln)
    file(strain_ids)

  output:
    file("FINAL*core.aln") optional true

  script:
    """
    for file in ${aln}; do
      NUM=\$(grep -c '>' \$file)
      STRAINS=\$(cat ${strain_ids} | wc -l)
      if [ \$NUM -eq \$STRAINS ]; then
        cp \${file} "\$(basename \${file} .aln)"_core.aln
        sed -r -i '/>/ s/_[^_]+\$//' "\$(basename \${file} .aln)"_core.aln
        sed -r -i '/^[^>]/ s/-/X/g' "\$(basename \${file} .aln)"_core.aln
      fi
    done

    for file in *_core.aln; do
      cd-hit -i \$file -o "\$file"_TMP -c 1.0
      RECORDS=\$(grep -c ">" "\$file"_TMP)
      if [ \$RECORDS -ne 1 ]; then
        cp \$file FINAL_\$file
        sed -r -i '/^[^>]/ s/X/-/g' FINAL_\$file
      fi      
    done

    """
}

