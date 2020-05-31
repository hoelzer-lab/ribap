#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
Nextflow -- RIBAP
Author: hoelzer.martin@gmail.com
        kevin.lamkiewicz@uni-jena.de
*/

/************************** 
* META & HELP MESSAGES 
**************************/

/* 
Comment section: First part is a terminal print for additional user information,
followed by some help statements (e.g. missing input) Second part is file
channel input. This allows via --list to alter the input of --nano & --illumina
to add csv instead. name,path   or name,pathR1,pathR2 in case of illumina 
*/

// terminal prints
if (params.help) { exit 0, helpMSG() }

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Workdir location:"
println "  $workflow.workDir\u001B[0m"
println " "
if (workflow.profile == 'standard') {
println "\033[2mCPUs to use: $params.cores"
println "Output dir name: $params.output\u001B[0m"
println " "}

if( !nextflow.version.matches('20.01+') ) {
    println "This workflow requires Nextflow version 20.01 or greater -- You are running version $nextflow.version"
    exit 1
}

if (params.profile) { exit 1, "--profile is WRONG use -profile" }
if (params.fasta == '' ) { exit 1, "input missing, use [--fasta]"}

/************************** 
* INPUT CHANNELS 
**************************/

// genome fasta input & --list support
if (params.fasta && params.list) { fasta_input_ch = Channel
  .fromPath( params.fasta, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
  .view() }
  else if (params.fasta) { fasta_input_ch = Channel
    .fromPath( params.fasta, checkIfExists: true)
    .map { file -> tuple(file.simpleName, file) }
}

/************************** 
* MODULES
**************************/

/* Comment section: */

include rename from './modules/rename' 
include prokka from './modules/prokka' 
include strain_ids from './modules/strain_ids' 
include roary from './modules/roary' 
include {mmseqs2; mmseq2tsv} from './modules/mmseqs2'
include ilp_build from './modules/ilp_build'
include ilp_solve from './modules/ilp_solve' 
include combine_roary_ilp from './modules/combine_roary_ilp'
include prepare_msa from './modules/prepare_msa' 
include mafft from './modules/mafft' 
include fasttree from './modules/fasttree'
include nw_display from './modules/nw_display' 
include combine_msa from './modules/combine_msa'
include generate_html from './modules/generate_html'
include generate_upsetr_input from './modules/generate_upsetr_input' 
include upsetr from './modules/upsetr' 
if (params.sets) {include upsetr_subset from './modules/upsetr'}
if (params.tree) {include raxml from './modules/raxml'}


/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {

  prokka(rename(fasta_input_ch))
  
  gff_ch = prokka.out[0]
  faa_ch = prokka.out[1].collect()

  strain_ids(prokka.out[0].collect())

  identity_ch = Channel.from(60, 70, 80, 90, 95)
  roary_run_ch = identity_ch.combine(gff_ch).groupTuple()
  roary(roary_run_ch)

  mmseqs2(faa_ch)

  ilp_solve(
    ilp_build(
      mmseq2tsv(mmseqs2.out, strain_ids.out).flatten()
    )
  )

  

  //copy all *sol and *simple into a solved folder for ilp_solve
  combine_ch = identity_ch
    .join(roary.out)
    .concat(strain_ids.out)
    .join(identity_ch
      .combine(strain_ids.out))
    .join(identity_ch
      .combine(gff_ch).groupTuple())

  combine_roary_ilp(combine_ch, ilp_solve.out[0].flatten().toList()) 



  // select only the 95 combined output file
  identity_ch = Channel.from(95)
  prepare_msa(identity_ch.join(combine_roary_ilp.out[0]), prokka.out[1].map { id, faa -> faa}.collect())

  // 50 alignments will be processed one after the other
  nw_display(
    fasttree(
      mafft(  
        prepare_msa.out.flatten().buffer(size: 50, remainder: true)
      )
    )
  )

  combine_msa(mafft.out.collect(), strain_ids.out)

  build_html_ch = identity_ch.join(combine_roary_ilp.out[0])
  generate_html(build_html_ch, roary.out.collect(), combine_roary_ilp.out[1].collect(), nw_display.out.collect())

  generate_upsetr_input(identity_ch.join(combine_roary_ilp.out[0]), strain_ids.out)
  upsetr(generate_upsetr_input.out[1])
  if (params.sets) {upsetr_subset(generate_upsetr_input.out[1])}

  if (params.tree) {raxml(combine_msa.out)}

}



/**************************  
* --help
**************************/
def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________
    
    RIBAP - Roary ILP Bacterial Annotation Pipeline

    Annotate your protein sequences with Prokka and determine a pan genome with Roary.
    This genome is refined with the usage of ILPs that solve the best matching for each pairwise
    strain blastp comparison.
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run ribap.nf --fasta '../strains/*.fasta' 

    ${c_yellow}Input:${c_reset}
    ${c_green} --fasta ${c_reset}           '*.fasta'         -> one strain per file
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset}            

    ${c_yellow}Params:${c_reset}
    --tmlim             Time limit for ILP solve [default: $params.tmlim]
    --gcode             Genetic code for Prokka annotation [default: $params.gcode]
    --tree              build tree based on the core genome? 
                        Sure thing, We will use RAxML for this. 
                        Be aware, this will take a lot of time. [default: $params.tree]

    ${c_yellow}UpSet plot:${c_reset}
    --sets              FASTA simpleNames for genomes that should be 
                        used in the UpSet plotting. Needed format:
                        "\\"Cav\\",\\"Cab\\",\\"Cga\\",\\"Ctr\\"" [default: $params.sets]
                        ${c_dim}(sorry, this will be simplified someday)${c_reset}
    --heigth            Height of the plot [default: $params.heigth]
    --width             Width of the plot [default: $params.width]

    ${c_yellow}Compute options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --memory            max memory for local use [default: $params.memory]
    --output            name of the result folder [default: $params.output]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    ${c_yellow}LSF computing:${c_reset}
    For execution of the workflow on a HPC with LSF adjust the following parameters:
    --databases         defines the path where databases are stored [default: $params.cloudDatabase]
    --workdir           defines the path where nextflow writes tmp files [default: $params.workdir]
    --cachedir          defines the path where images (singularity) are cached [default: $params.cachedir] 


    ${c_yellow}Profile:${c_reset}
    -profile                 standard (local, pure docker) [default]
                             conda (mixes conda and docker)
                             lsf (HPC w/ LSF, singularity/docker)
                             ebi (HPC w/ LSF, singularity/docker, preconfigured for the EBI cluster)
                             ${c_reset}
    """.stripIndent()
}

  