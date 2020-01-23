#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
Nextflow -- Analysis Pipeline
Author: hoelzer.martin@gmail.com
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

if (params.profile) { exit 1, "--profile is WRONG use -profile" }
if (params.db == '' &&  params.names == '' ) { exit 1, "input missing, use [--db] or [--names]"}

/************************** 
* INPUT CHANNELS 
**************************/

db = file(params.db)

names = file(params.names)
    .readLines()
    .each{it}


/************************** 
* MODULES
**************************/

/* Comment section: */

include sketch from './modules/mash' params(ksize: params.ksize)
include screen from './modules/mash' params(output: params.output)


/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {
  reads = Channel.fromSRA(names, apiKey: params.key)

  // Create a MinHash signature/ sketch of the reference collection
  sketch(db)

  // Search k-mers of read set against this signature
  screen(sketch.out, reads)
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
    
    Workflow: Template
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run wf_template --nano '*/*.fastq' 

    ${c_yellow}Input:${c_reset}
    ${c_green} --db ${c_reset}            fasta [default: $params.db]
    ${c_green} --names ${c_reset}         txt [default: $params.names]

    ${c_yellow}Options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --memory            max memory for local use [default: $params.memory]
    --output            name of the result folder [default: $params.output]

    ${c_yellow}Parameters:${c_reset}
    --key             NCBI ket [default: $params.key]
    --ksize           mash kmer size [default: $params.ksize]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    ${c_yellow}LSF computing:${c_reset}
    For execution of the workflow on a HPC with LSF adjust the following parameters:
    --databases         defines the path where databases are stored [default: $params.cloudDatabase]
    --workdir           defines the path where nextflow writes tmp files [default: $params.workdir]
    --cachedir          defines the path where images (singularity) are cached [default: $params.cachedir] 


    Profile:
    -profile                 standard (local, pure docker) [default]
                             conda (mixes conda and docker)
                             lsf (HPC w/ LSF, singularity/docker)
                             ebi (HPC w/ LSF, singularity/docker)
                             ${c_reset}
    """.stripIndent()
}

  