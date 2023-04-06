/*Comment section: */

process ilp_refinement {
  label 'python3' 
  publishDir "${params.output}/04-ilp", pattern: "*.ilp.simple" 

  input: 
    file(tsv)

  output:
  path("*.ilp.simple")


  script:
    def KEEP = params.keepILPs ? "--keep" : ""
    """
    derive_ilp_solutions.py --tmlim ${params.tmlim} --max --indel ${tsv} ${KEEP}
    """
}

