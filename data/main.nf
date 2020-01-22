nextflow.preview.dsl = 2


// Gimme all you got
threads = 8

// The fasta reference genomes we want to query reads against, e.g. 
// VIPR CoV genomes or a single genome of the Wuhan 2019-nCoV strain
db = file("47907236622-GenomicFastaResults.fasta")

// Where to store results
outdir = "foobar"

// We want a list of SRA IDs like names = ["SRR7548061", "SRR7548062"]
// cut -f1 -d, info.csv | tail -n +2 > info.txt
names = file('info.txt')
    .readLines()
    .each{it}

// The API key so NCBI SRA does not get upset when we want to download a lot
key = "ed38c182ff2cdca964f0c766e220c02ec608"


workflow {
    // https://www.nextflow.io/docs/edge/channel.html#fromsra
    // https://www.nextflow.io/blog/2019/release-19.03.0-edge.html
    // https://github.com/nextflow-io/nextflow/issues/1196
    // https://github.com/nanozoo/wf_metagenomics/issues/15
    // println(names)
    reads = Channel.fromSRA(names, apiKey: key, max: 10)
    // reads.view()
    // [ERR3569479, [.../ERR3569479_1.fastq.gz, .../ERR3569479_2.fastq.gz]]
    
    // Create a MinHash signature/ sketch of the reference collection
    sketch(db)

    // Search k-mers of read set against this signature
    screen(sketch.out, reads)
}



process sketch {
    input:
    file(db)

    output:
    file("db.msh")

    """
    mash sketch -i -s 1000 -k 15 -S 42 -o db.msh ${db}
    """
}


process screen {
    publishDir "${outdir}", mode: "copy", pattern: "${uid}.screen.csv"
    
    input:
    file(sketch)
    tuple val(uid), file(reads)

    output:
    // file "${uid}.screen.csv"
    file("${uid}.screen.csv")

    """
    echo ${reads}
    mash screen -w -p ${threads} -i 0 ${sketch} ${reads} > ${uid}.screen.csv
    """
    // Streaming from 2 inputs ...
    // "reads" is a list of all fastq files nf finds for each SRA ID,
    // and all will be fed into one call of mash -- great stuff.
}