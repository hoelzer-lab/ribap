/*Comment section: */

process ilp_refinement {
  label 'python3' 
  publishDir "${params.output}/05-ilp", pattern: "*.ilp.simple" 

  input: 
    file(tsv)

  output:
  path("*.ilp.simple")


  script:
    def KEEP = params.keepILPs ? "--keep" : ""
    """
    derive_ilp_solutions.py -p ${task.cpus} --tmlim ${params.tmlim} --max --indel ${tsv} ${KEEP}
    """
}

