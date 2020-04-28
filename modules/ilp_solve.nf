/*Comment section: */

process ilp_solve {
  label 'glpk'
//  publishDir "${params.output}/ilp/solved", mode: 'copy', pattern: "*sol" 
  publishDir "${params.output}/ilp/solved", mode: 'copy', pattern: "simple*" 

  input: 
    tuple val(name), path(ilp)

  output:
//    tuple file("solved/*.sol"), file("simple/*.simple")
      path("simple*", type: 'dir')

  script:
    """
    mkdir solved
    for ILP in ilp/*.ilp; do 
      BN=\$(basename \$ILP .ilp)
      glpsol --lp \$ILP --mipgap 0.01 --pcost --cuts --memlim 16834 --tmlim ${params.tmlim} -o solved/\$BN.sol  
      sed -E -i '/x_A[^[:space:]]+\$/ N;s/\\n//g' "solved/\$BN.sol"
    done

    TMP=\$(basename \$PWD)
    mkdir simple_"\$TMP"
    for SOL in solved/*.sol; do
        awk '\$2 ~ /x_A.*_/ && \$4 == 1 {print}' "\$SOL"
    done > simple_"\$TMP"/"${name}".ilp.simple
    """
}
