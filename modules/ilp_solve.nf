/*Comment section: */

process ilp_solve {
  label 'glpk'
  publishDir "${params.output}/ilp/", mode: 'copy', pattern: "solved/*.sol" 
  publishDir "${params.output}/ilp/simple", mode: 'copy', pattern: "simple*" 

// there is a problem with multiple use of variable 'x_A1h_A1t', see issue #11
// it seems that the problem only occurs for longer tmlim because the 33er set was running 
// through with 240s but the error occured for 7200s
  errorStrategy{task.exitStatus=101 ?'ignore':'terminate'}

  input: 
    tuple val(name), path(ilp)

  output:
    path("simple*", type: 'dir')
    path("solved/*.sol")

  script:
    """
    # can we use parallel inside a docker? seems so
    mkdir solved
    cp ilp/*.ilp .
    ls *.ilp | parallel -j "${task.cpus}" -I% --max-args 1 glpsol --lp % --mipgap 0.01 --pcost --cuts --memlim 16834 --tmlim ${params.tmlim} -o solved/%.sol
    rm *.ilp
    for SOL in solved/*.sol; do
      sed -E -i '/x_A[^[:space:]]+\$/ N;s/\\n//g' \$SOL
    done

    # this is the succesive version w/o parallel
#    mkdir solved
#    for ILP in ilp/*.ilp; do 
#      BN=\$(basename \$ILP .ilp)
#      glpsol --lp \$ILP --mipgap 0.01 --pcost --cuts --memlim 16834 --tmlim ${params.tmlim} -o solved/\$BN.sol  
#      sed -E -i '/x_A[^[:space:]]+\$/ N;s/\\n//g' "solved/\$BN.sol"
#    done

    TMP=\$(basename \$PWD)
    mkdir simple_"\$TMP"
    for SOL in solved/*.sol; do
        awk '\$2 ~ /x_A.*_/ && \$4 == 1 {print}' "\$SOL"
    done > simple_"\$TMP"/"${name}".ilp.simple
    """
}
