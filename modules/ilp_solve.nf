/*Comment section: */

process ilp_solve {
  label 'glpk'
  //publishDir "${params.output}/ilp/", pattern: "solved/*.sol" 
  publishDir "${params.output}/05-ilp", pattern: "simple*/*.simple" 

//  there was a problem with multiple use of variable 'x_A1h_A1t', see issue #11
//  errorStrategy{task.exitStatus=101 ?'ignore':'terminate'}

  input: 
    tuple val(name), path(ilp), path(files)

  output:
    path("simple*", type: 'dir')
    //path("solved/*.sol")
    path("simple*/*.simple")

  script:

    def delete = "${params.deleteILPs}"

    """
    # can we use parallel inside a docker? seems so
    mkdir solved
    ls *.ilp | parallel -j "${task.cpus}" --max-args 1 glpsol --lp {} --mipgap 0.01 --pcost --cuts --memlim 16834 --tmlim ${params.tmlim} -o solved/{}.sol \
                        && if [[ ${delete} == true ]] ; then rm {}; fi

    TMP=\$(basename ilp_*)
    mkdir simple_"\$TMP"

    if [[ ${delete} == true ]] ; then 
      rm -r \$(realpath ${ilp})
    fi

    for SOL in solved/*.sol; do
      sed -i '/x_A[^[:space:]]+\$/ N;s/\\n//g' \$SOL | awk '\$2 ~ /x_A.*_/ && \$4 == 1 {print}' 
      if [[ ${delete} == true ]] ; then 
        rm \$SOL
      fi
    done > simple_"\$TMP"/"${name}".ilp.simple
    
    """
}
