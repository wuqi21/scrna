process PARSE_PROTOCOL {
    tag "$meta.id"
    label 'process_single'

    conda 'bioconda::pyfastx=2.1.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyfastx:2.1.0--py39h3d4b85c_0' :
        'biocontainers/pyfastx:2.1.0--py39h3d4b85c_0' }"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), path(reads)
    path json_file
    val protocol

    output:
    path "*starsolo_cb_umi_args.txt", emit: starsolo_cb_umi_args
    path "*parsed_protocol.txt", emit: parsed_protocol
    path "*whitelist.txt", emit: whitelist

    when:
    task.ext.when == null || task.ext.when

    script:
    def new_protocol = protocol == 'new' ? "" : ""
    """
    parse_protocol.py \\
        --sample ${meta.id} \\
        --R1_fastq ${reads[0]} \\
        --json_file ${json_file} \\
        --protocol ${protocol} \\
        ${new_protocol}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_protocol.py: \$(parse_protocol.py --version)
    END_VERSIONS
    """
}
