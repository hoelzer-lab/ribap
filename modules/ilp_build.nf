/*Comment section: */

process ilp_build {
  label 'python3'
  publishDir "${params.output}/ilp", mode: 'copy', pattern: "ilp" 

  input: 
    file(tsv)

  output:
    tuple env(BN), path('ilp', type: 'dir')

  script:
    """
    #parallel -j "${task.cpus}" 'ILP.py --max --indel {} > '"\$PWD"'/{/.}.ilp 2>/dev/null' ::: "\$PWD"/*tsv 2>/dev/null
    #ILP.py --max --indel ${tsv} > \$(basename ${tsv} .tsv).ilp
    
    BN=\$(basename ${tsv} .tsv)

    mkdir ilp
    ILP.py --max --indel ${tsv}

    echo \$BN
    """
}

