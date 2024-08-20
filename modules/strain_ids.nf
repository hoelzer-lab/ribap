/*Comment section: */

process strain_ids {
  label 'basics'
  publishDir "${params.output}/", mode: 'copy', pattern: "strain_ids.txt" 

  input: 
    file(gff)

  output:
    file("strain_ids.txt")

  script:
    """
    for ANNO in *.gff*; do
      BN=\$(basename "\$ANNO" .gff)
      echo \$(grep -vm1 '^#' \$ANNO | awk '{print \$9}' | cut -d'=' -f2 | cut -d'_' -f1),"\$BN"
    done > strain_ids.txt
    """
}

