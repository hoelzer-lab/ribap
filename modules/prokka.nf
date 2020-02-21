/*Comment section: */

process prokka {
  label 'prokka'
  publishDir "${params.output}/prokka", mode: 'copy', pattern: "${name}.gff" 
  publishDir "${params.output}/prokka", mode: 'copy', pattern: "${name}.faa" 

  input: 
    tuple val(name), file(fasta)

  output:
    file("${name}.gff")
    tuple val(name), file("${name}.faa")

  script:
    """
      prokka --gcode ${params.gcode} --cpus ${task.cpus} --outdir output --prefix annotation ${fasta}
      mv output/annotation.faa ${name}.faa
      mv output/annotation.gff ${name}.gff
      mv output/annotation.gbk ${name}.gbk    
    """
}

