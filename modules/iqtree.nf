/*Comment section: */

process iqtree {
  label 'iqtree'
  publishDir "${params.output}/iqtree", mode: 'copy', pattern: "*modeltest*" 

  input: 
    file(aln)
    file(nexus)

  output:
    file("*modeltest*")

  script:
    """
    iqtree -p ${nexus} -T AUTO --threads-max ${task.cpus} -nt 4 -m TEST -pre "\$(basename ${nexus} .nex)"-modeltest
    """
}