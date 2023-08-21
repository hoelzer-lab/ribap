/*Comment section: */

process filter_alignment {
  label 'filter_alignment'
  publishDir "${params.output}/07-mafft", pattern: "FINAL*core.aln" 

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

      # calculate cutoff for how many species are needed to define a core gene
      # default is 1.0, but can be lowered by the user
      t=\$(echo \$STRAINS*$params.core_perc | bc)
      # round, everything equal or below .5 will be rounded down, otherwise up
      # printf "%.0f" 26.4 == 26
      # printf "%.0f" 26.52 == 27
      STRAINS_CUTOFF=\$(printf "%.0f" \$t)

      if [ \$NUM -ge \$STRAINS_CUTOFF ]; then
        cp \${file} "\$(basename \${file} .aln)"_core.aln
        sed -r -i '/>/ s/_[^_]+\$//' "\$(basename \${file} .aln)"_core.aln
        sed -r -i '/^[^>]/ s/-/X/g' "\$(basename \${file} .aln)"_core.aln
      fi
    done

    # it can happen that there is no MSA with all input species!
    if ls ./ | grep '_core.aln'; then
      for file in *_core.aln; do
        cd-hit -i \$file -o "\$file"_TMP -c 1.0
        RECORDS=\$(grep -c ">" "\$file"_TMP)
        if [ \$RECORDS -ne 1 ]; then
          cp \$file FINAL_\$file
          sed -r -i '/^[^>]/ s/X/-/g' FINAL_\$file
        fi      
      done
    fi
    """
}

