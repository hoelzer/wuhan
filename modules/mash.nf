/*Comment section: */


process sketch {
    label 'mash'  

    input:
    file(db)

    output:
    file("db.msh")

    """
    mash sketch -p ${task.cpus} -i -s 1000 -k 15 -S 42 -o db.msh ${db}
    """
}


process screen {
    label 'mash'  
    publishDir "${params.output}/", mode: 'copy', pattern: "${uid}.screen.csv"
    
    input:
    file(sketch)
    tuple val(uid), file(reads)

    output:
    // file "${uid}.screen.csv"
    file("${uid}.screen.csv")

    """
    echo ${reads}
    mash screen -w -p ${task.cpus} -i 0 ${sketch} ${reads} > ${uid}.screen.csv
    """
    // Streaming from 2 inputs ...
    // "reads" is a list of all fastq files nf finds for each SRA ID,
    // and all will be fed into one call of mash -- great stuff.
}


