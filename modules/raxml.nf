/*Comment section: */

process raxml {
  label 'raxml'
  publishDir "${params.output}/raxml", mode: 'copy', pattern: "*.aa" 

  input: 
    file(aln)

  output:
    file("coreGenome_mafft.raxml.support")

  script:
    """
    #raxmlHPC-PTHREADS-SSE3 -T ${task.cpus} -f a -x 1234 -p 1234 -s coreGenome_mafft.aln -n aa -m PROTGAMMAWAG -N 100 
    raxml-ng --all --threads ${task.cpus} --msa coreGenome_mafft.aln --prefix coreGenome_mafft --msa-format FASTA --model PROTGTR+G --bs-trees 100 
    """
}

