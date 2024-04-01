process STARSOLO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'biocontainers/star:2.7.11b--h43eeafb_0' }"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), path(reads)
    path index
    path assets_dir
    val starsolo_cmd

    output:
    tuple val(meta), path("${meta.id}.matrix/")       , emit: matrix
    tuple val(meta), path('*d.out.bam')       , emit: bam
    tuple val(meta), path('*.Solo.out/*')       , emit: counts
    path('*Log.final.out')   , emit: log_final
    path "*.Solo.out/GeneFull_Ex50pAS/${meta.id}.Summary.csv" , emit: summary
    path  "versions.yml"                      , emit: versions

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*out.mate')               , optional:true, emit: unmap
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"

    // do not indent
"""
${starsolo_cmd}

if [ -d ${prefix}.Solo.out ]; then
    # Backslashes still need to be escaped (https://github.com/nextflow-io/nextflow/issues/67)
    find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip -f {} \\;
fi

mv ${prefix}.Solo.out/GeneFull_Ex50pAS/Summary.csv ${prefix}.Solo.out/GeneFull_Ex50pAS/${prefix}.Summary.csv
mkdir ${prefix}.matrix
mv ${prefix}.Solo.out/GeneFull_Ex50pAS/{raw,filtered} ./${prefix}.matrix/

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    star: \$(STAR --version | sed -e "s/STAR_//g")
END_VERSIONS
"""
}
