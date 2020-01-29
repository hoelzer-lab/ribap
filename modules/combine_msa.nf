/*Comment section: */

process combine_msa {
  label 'pyhton3'
  publishDir "${params.output}/msa", mode: 'copy', pattern: "coreGenome_mafft.aln" 

  input: 
    file(aln)
    file(strain_ids)

  output:
    file("coreGenome_mafft.aln")

  script:
    """
    concat_coreMSA.py . \$(wc -l ${strain_ids} | awk '{print \$1}')
    """
}

