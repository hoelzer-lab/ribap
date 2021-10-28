/*Comment section: */

process prokka {
  label 'prokka'
  publishDir "${params.output}/prokka", mode: 'copy', pattern: "${name}" 

  input: 
    tuple val(name), file(fasta)

  output:
    file("${name}/${name}.gff")
    tuple val(name), file("${name}/${name}.faa")
    path("${name}", type: 'dir')

  script:
    """
      prokka --gcode ${params.gcode} --cpus ${task.cpus} --outdir ${name} --prefix ${name} ${fasta}
    """
}

