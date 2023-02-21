#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
channel input. This allows via --list to alter the input of --fasta
to add csv instead. name,path  
*/

// terminal prints
if (params.help) { exit 0, helpMSG() }

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Workdir location (intermediate files):"
println "  $workflow.workDir"
println "Output dir name:"
println "  $params.output\u001B[0m"
println " "
if (workflow.profile == 'standard' || workflow.profile == 'local') {
println "\033[2mCPUs to use: $params.cores"}
if ( workflow.profile.contains('singularity') ) {
    println "\033[2mSingularity cache directory:"
    println "  $params.singularityCacheDir\u001B[0m"
}
if ( workflow.profile.contains('conda') ) { 
    println "\033[2mConda cache directory:"
    println "  $params.condaCacheDir\u001B[0m"
}
if (params.profile) { exit 1, "--profile is WRONG use -profile" }
if (params.fasta == '' ) { exit 1, "input missing, use [--fasta]"}

if ( !workflow.revision ) { 
    println ""
    println "\033[0;33mWARNING: Not a stable execution. Please use -r for full reproducibility. Use nextflow info hoelzer-lab/ribap to list release versions.\033[0m"
    println "\033[0;33mIf you cloned the github repository, make sure to use the latest stable release or be aware of potential bugs.\033[0m"
}

def folder = new File(params.output)
if ( folder.exists() ) { 
    println ""
    println "\033[0;33mWARNING: Output folder already exists. Results might be overwritten! You can adjust the output folder via [--output]\033[0m"
}

if ( workflow.profile.contains('singularity') ) {
    println ""
    println "\033[0;33mWARNING: Singularity image building sometimes fails!"
    println "Multiple resumes (-resume) and --max_cores 1 --cores 1 for local execution might help.\033[0m"
}

if ( params.bootstrap < 1000 ) { exit 1, "--bootstrap needs to be >=1000 (IQ-TREE -bb parameter requirement for ultra-fast bootstraping)"}

if ( params.deleteILPs ) { 
    println ""
    println "\033[0;33mWARNING: ILPs will be deleted on the fly which means you can not -resume them (but other steps of the pipeline).\033[0m"
}

/************************** 
* INPUT CHANNELS 
**************************/

// genome fasta input & --list support
if (params.fasta && params.list) { fasta_input_ch = Channel
  .fromPath( params.fasta, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
  //.view() 
  }
  else if (params.fasta) { fasta_input_ch = Channel
    .fromPath( params.fasta, checkIfExists: true)
    .map { file -> tuple(file.baseName, file) }
}

// reference gbk file to improve annotation & --list support 
// list support expects something like:
//
// genome1,/path/to/ref1.gbk
// genome2,/path/to/ref1.gbk
//
// where the first column matches the basenames of the input genome FASTAs
if (params.reference && params.list) { reference_input_ch = Channel
  .fromPath( params.reference, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
  //.view() 
  }
  else if (params.reference) { reference_input_ch = Channel
    .fromPath( params.reference, checkIfExists: true)
    .map { file -> file }
} else {
  reference_input_ch = Channel.value('null')
}


/************************** 
* MODULES
**************************/

/* Comment section: */

include { rename } from './modules/rename' 
include { prokka } from './modules/prokka' 
include { strain_ids } from './modules/strain_ids' 
include { roary } from './modules/roary' 
include { mmseqs2; mmseqs2tsv } from './modules/mmseqs2'
include { ilp_build } from './modules/ilp_build'
include { ilp_solve } from './modules/ilp_solve' 
include { combine_roary_ilp } from './modules/combine_roary_ilp'
include { prepare_msa } from './modules/prepare_msa' 
include { mafft } from './modules/mafft' 
include { fasttree } from './modules/fasttree'
include { nw_display } from './modules/nw_display' 
// include { combine_msa } from './modules/combine_msa'
include { generate_html } from './modules/generate_html'
include { generate_upsetr_input } from './modules/generate_upsetr_input' 
include { upsetr } from './modules/upsetr' 
if (params.sets) {include { upsetr_subset } from './modules/upsetr'}

if (params.tree) {
  include { filter_alignment } from './modules/filter_alignment'
  include { nexus_core } from './modules/nexus_core'
  include { iqtree } from './modules/iqtree'
  }


/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {

  renamed_fasta_ch = rename(fasta_input_ch)

  if (params.reference && params.list) {
    prokka_input_ch = renamed_fasta_ch.join(reference_input_ch, remainder: true).map { id, id_renamed, fasta, gbk -> [id_renamed, fasta, gbk]}
  } else {
    // this will either produce a channel w/ [sample_RENAMED, fasta_RENAMED, reference_gbk] OR [sample_RENAMED, fasta_RENAMED, null] 
    prokka_input_ch = renamed_fasta_ch.combine(reference_input_ch).map { id, id_renamed, fasta, gbk -> [id_renamed, fasta, gbk]}
  }

  prokka(prokka_input_ch)

  gff_ch = prokka.out[0]
  faa_ch = prokka.out[1].collect()

  strain_ids(prokka.out[0].collect())

  // identity_ch = Channel.from(60, 70, 80, 90, 95)
  // roary_run_ch = identity_ch.combine(gff_ch).groupTuple()
  // roary(roary_run_ch)

  mmseqs2(faa_ch)

  // mmseqs2tsv(mmseqs2.out[0], strain_ids.out)
   ilp_solve(
    ilp_build(
      mmseqs2tsv(mmseqs2.out[0], strain_ids.out).flatten()
    )
  )
  


//   // select only the 95 combined output file
//   identity_ch = Channel.from(95)
//   //copy all *sol and *simple into a solved folder for ilp_solve
//   combine_ch = identity_ch
//     .join(roary.out)
//     .concat(strain_ids.out)
//     .join(identity_ch
//       .combine(strain_ids.out))
//     .join(identity_ch
//       .combine(gff_ch).groupTuple())

//   combine_roary_ilp(combine_ch, ilp_solve.out[0].flatten().toList()) 



//   // // select only the 95 combined output file
//   // identity_ch = Channel.from(95)
//   prepare_msa(identity_ch.join(combine_roary_ilp.out[0]), prokka.out[1].map { id, faa -> faa}.collect())

//   // 50 alignments will be processed one after the other
//   nw_display(
//     fasttree(
//       mafft(  
//         prepare_msa.out.flatten().buffer(size: 50, remainder: true)
//       )
//     )
//   )

//   //combine_msa(mafft.out.collect(), strain_ids.out)
//   build_html_ch = identity_ch.join(combine_roary_ilp.out[0])
//   generate_html(build_html_ch, roary.out.collect(), combine_roary_ilp.out[1].collect(), nw_display.out.collect())

//   generate_upsetr_input(identity_ch.join(combine_roary_ilp.out[0]), strain_ids.out)
//   upsetr(generate_upsetr_input.out[1])
//   if (params.sets) {upsetr_subset(generate_upsetr_input.out[1])}

//   //if (params.tree) {raxml(combine_msa.out)}
//   if (params.tree) {
//     filter_alignment(mafft.out.collect(), strain_ids.out)
//     nexus_core(filter_alignment.out.collect())
//     iqtree(filter_alignment.out.collect(), nexus_core.out)
//   }

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

    Annotate your bacterial genome sequences with Prokka and determine a pangenome with Roary.
    The Roary gene clusters are further refined with the usage of ILPs that solve the best matching 
    for each pairwise strain MMSeqs2 comparison.
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run ribap.nf --fasta '../strains/*.fasta' 

    ${c_yellow}Input:${c_reset}
    ${c_green} --fasta ${c_reset}           '*.fasta'         -> one strain per file
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset}            

    ${c_yellow}Params:${c_reset}
    --tmlim             Time limit for ILP solve [default: $params.tmlim]
    --gcode             Genetic code for Prokka annotation [default: $params.gcode]
    --reference         A reference genbank (gbk, gb) file to guide functional annotation via Prokka.
                        Attention: when directly provided without the --list parameter, all input genomes 
                        will be functionally annotated using the same reference. To use different reference files
                        or to exclude certain genomes from reference-based annotation use the --list option. [defaut: $params.reference]
    --tree              build tree based on the core genome? 
                        Sure thing, We will use RAxML for this. 
                        Be aware, this will take a lot of time. [default: $params.tree]
    --bootstrap         Bootstrap value for tree building (increases time!). Must be >=1000 for IQ-TREE ultra-fast bootstraps [default: $params.bootstrap] 

    ${c_yellow}UpSet plot:${c_reset}
    --sets              FASTA simpleNames for genomes that should be 
                        used in the UpSet plotting. Needed format:
                        "\\"Cav\\",\\"Cab\\",\\"Cga\\",\\"Ctr\\"" [default: $params.sets]
                        ${c_dim}(sorry, this will be simplified someday)${c_reset}
    --heigth            Height of the plot [default: $params.heigth]
    --width             Width of the plot [default: $params.width]

    ${c_yellow}Compute options:${c_reset}
    --cores             max cores used per process for local use [default: $params.cores]
    --max_cores         max cores used on the machine for local use [default: $params.max_cores]
    --memory            max memory for local use [default: $params.memory]
    --output            name of the result folder [default: $params.output]
    --deleteILPs        the ILPs can take a lot (!) of space. Use this flag to delete them on the fly.
                        Attention: you can not -resume already calculated ILPs when deleting them. [default: $params.deleteILPs]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    ${c_yellow}Caching:${c_reset}
    --condaCacheDir          Location for storing the conda environments [default: $params.condaCacheDir]
    --singularityCacheDir    Location for storing the singularity images [default: $params.singularityCacheDir]
    -w                	     Working directory for all intermediate results [default: $workflow.workDir]

    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
    
    ${c_green}Executer${c_reset} (choose one):
      local
      slurm
      lsf
    
    ${c_blue}Engines${c_reset} (choose one):
      conda
      docker
      singularity
    
    Per default: -profile local,conda is executed. 
    """.stripIndent()
}

  