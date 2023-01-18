/*Comment section: */

process raxml {
  label 'raxml'
  publishDir "${params.output}/raxml", mode: 'copy', pattern: "*.support" 

  input: 
    file(aln)

  output:
    file("coreGenome_mafft.raxml.support")

  script:
    """
    #raxmlHPC-PTHREADS-SSE3 -T ${task.cpus} -f a -x 1234 -p 1234 -s ${aln} -n aa -m PROTGAMMAWAG -N ${params.bootstrap}
    raxml-ng --all --threads auto{${task.cpus}} --msa ${aln} --prefix coreGenome_mafft --msa-format FASTA --model PROTGTR+G --bs-trees ${params.bootstrap} 
    """
}

// --threads auto{4} defined based on https://github.com/amkozlov/raxml-ng/wiki/Parallelization#core-oversubscription and bc/ error on SLURM HPC