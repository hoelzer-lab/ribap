/*Comment section: */

process ilp_build {
  label 'python3'
  //publishDir "${params.output}/ilp/sub_ilps", pattern: "ilp*/*.ilp" 

  input: 
    file(tsv)

  output:
  // tuple path("ilp_*", type: 'dir'), path("ilp_*/*.ilp")
  path("ilp_*", type: 'dir')
  path("ilp_*/*.ilp")


  script:
    """
    #parallel -j "${task.cpus}" 'ILP.py --max --indel {} > '"\$PWD"'/{/.}.ilp 2>/dev/null' ::: "\$PWD"/*tsv 2>/dev/null
    
    #BN=\$(basename ${tsv} .pkl)

    #mkdir ilp
    derive_ilp_solutions.py -p ${params.cores} --tmlim ${params.tmlim} --max --indel ${tsv}
    #for BN in \$(for file in ilp/*ilp; do echo \${file%_*}; done | sort | uniq); do
    #  BN=\$(basename \$BN)
    #  mkdir ilp_"\${BN}"
    #  mv ilp/"\${BN}"*.ilp ilp_"\${BN}" 
    #done

    #rm -r ilp/

    #echo \$BN
    """
}

