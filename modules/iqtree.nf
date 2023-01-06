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
    iqtree -p ${nexus} -bb 1000 --threads-max ${task.cpus} -nt AUTO -m TEST -pre "\$(basename ${nexus} .nex)"-modeltest
    """
}