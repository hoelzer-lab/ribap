/*Comment section: */

process ilp_build {
  label 'python3' 
  publishDir "${params.output}/05-ilp", pattern: "*.ilp.simple" 

  input: 
    file(tsv)

  output:
  path("*.ilp.simple")


  script:
    """
    derive_ilp_solutions.py -p ${params.cores} --tmlim ${params.tmlim} --max --indel ${tsv}
    """
}

