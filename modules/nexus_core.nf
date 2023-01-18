/*Comment section: */

process nexus_core {
  label 'nexus'
  publishDir "${params.output}/12-nexus", mode: 'copy', pattern: "core_genome.nex" 

  input: 
    file(aln)

  output:
    file("core_genome.nex")

  script:
    """
    echo '#nexus' > core_genome.nex
    echo 'begin sets;' >> core_genome.nex
    i=1
    for file in ${aln}; do
      echo "charset part\$i = \$file: *;" >> core_genome.nex
      i="\$((i+1))"
    done
    echo 'end;' >> core_genome.nex

    """
}