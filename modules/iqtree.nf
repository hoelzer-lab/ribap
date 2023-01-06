/*Comment section: */

process iqtree {
  label 'iqtree'
  publishDir "${params.output}/iqtree", mode: 'copy', pattern: "*modeltest" 

  input: 
    file(aln)

  output:
    file("*modeltest")

  script:
    """
    iqtree -s ${aln} --mem ${task.memory} -T AUTO --threads-max ${task.cpus} -nt 4 -m TEST -pre "\$(basename \${aln})"-modeltest
    """
}