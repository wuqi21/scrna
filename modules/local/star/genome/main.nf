process STAR_GENOME {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.10b--h9ee0642_0':
        'biocontainers/star:2.7.10b--h9ee0642_0' }"

    input:
    path fasta
    path gtf

    output:
    path("star") , emit: index
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def memory      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def include_gtf = gtf ? "--sjdbGTFfile $gtf" : ''
    def fasta_sa = ( Math.log(fasta.size()) / Math.log(2) ) / 2 - 1
    def sa = Math.floor( Math.min(14, fasta_sa) )
    """
    mkdir star
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        $include_gtf \\
        --runThreadN $task.cpus \\
        --genomeSAindexNbases ${sa} \\
        $memory \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """

    stub:
    if (gtf) {
        """
        mkdir star
        touch star/Genome
        touch star/Log.out
        touch star/SA
        touch star/SAindex
        touch star/chrLength.txt
        touch star/chrName.txt
        touch star/chrNameLength.txt
        touch star/chrStart.txt
        touch star/exonGeTrInfo.tab
        touch star/exonInfo.tab
        touch star/geneInfo.tab
        touch star/genomeParameters.txt
        touch star/sjdbInfo.txt
        touch star/sjdbList.fromGTF.out.tab
        touch star/sjdbList.out.tab
        touch star/transcriptInfo.tab

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    } else {
        """
        mkdir star
        touch star/Genome
        touch star/Log.out
        touch star/SA
        touch star/SAindex
        touch star/chrLength.txt
        touch star/chrName.txt
        touch star/chrNameLength.txt
        touch star/chrStart.txt
        touch star/genomeParameters.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    }
}
