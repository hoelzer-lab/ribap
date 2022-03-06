/*Comment section: */

process prokka {
  label 'prokka'
  publishDir "${params.output}/prokka", mode: 'copy', pattern: "${name}" 

  input: 
    tuple val(name), file(fasta), file(reference)

  output:
    file("${name}/${name}.gff")
    tuple val(name), file("${name}/${name}.faa")
    path("${name}", type: 'dir')

  script:
    """
    if [[ ${reference} == 'null' ]]; then
      prokka --gcode ${params.gcode} --cpus ${task.cpus} --outdir ${name} --prefix ${name} ${fasta}
    else
      prokka --gcode ${params.gcode} --cpus ${task.cpus} --outdir ${name} --prefix ${name} --proteins ${reference} ${fasta}
    fi
    """
}
