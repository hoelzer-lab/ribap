/*Comment section: */

process iqtree {
  label 'iqtree'
  publishDir "${params.output}/11-iqtree", mode: 'copy', pattern: "*modeltest*" 

  input: 
    file(aln)
    file(nexus)

  output:
    file("*modeltest*")

  script:
    """
    iqtree -spp ${nexus} -bb ${params.bootstrap} --threads-max ${task.cpus} -nt AUTO -m TEST -pre "\$(basename ${nexus} .nex)"-modeltest
    """
}