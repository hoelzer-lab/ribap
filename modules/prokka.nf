/*Comment section: */

process prokka {
  label 'prokka'
  publishDir "${params.output}/01-prokka", mode: 'copy', pattern: "${name}" 

  input: 
    tuple val(name), file(fasta), file(reference)

  output:
    file("${name}/${name}.gff")
    tuple val(name), file("${name}/${name}.faa")
    path("${name}", type: 'dir')

  script:
    """
    # this might be dangerous if the user provides some empty file and THINKS he provided some decent GBK file: this will still work w/o letting the user know
    if [[ -s ${reference} && \$(cat ${reference}) != 'null' ]]; then
      prokka --gcode ${params.gcode} --cpus ${task.cpus} --outdir ${name} --prefix ${name} --proteins ${reference} ${fasta}
    else
      prokka --gcode ${params.gcode} --cpus ${task.cpus} --outdir ${name} --prefix ${name} ${fasta}
    fi
    """
}
