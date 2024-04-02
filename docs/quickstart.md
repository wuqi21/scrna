# Quick Start

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.
> After that, you can [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Run the pipeline with test data using

```
nextflow run singleron-RD/scrna -profile test,docker --outdir ./results
```

The test config file can be found [here](../conf/test.config). This is equivalent to the following command

```
nextflow run singleron-RD/scrna \
 -profile docker \
 --outdir ./results \
 --input 'https://raw.githubusercontent.com/singleron-RD/scrna_test_data/master/GEXSCOPE-V2/samplesheet.csv' \
 --fasta 'https://github.com/nf-core/test-datasets/raw/scrnaseq/reference/GRCm38.p6.genome.chr19.fa' \
 --gtf 'https://github.com/nf-core/test-datasets/raw/scrnaseq/reference/gencode.vM19.annotation.chr19.gtf' \
 --genome_name 'mouse_chr19' \
 --keep_attributes 'gene_type=protein_coding,lncRNA,antisense' \
 --max_cpus 2 \
 --max_memory '6.GB' \
 --max_time '6.h'
```

If you prefer a web-based graphical interface or an interactive command-line wizard tool to enter the pipeline parameters, you can use [nf-core launch](https://oldsite.nf-co.re/tools/#launch-a-pipeline):

```
pip install nf-core
nf-core launch singleron-RD/scrna
```

Since indexing is an expensive process in time and resources you should ensure that it is only done once, by retaining the indices generated from each batch of reference files. When running data from the same genome later, you can provide the star genome index path to skip the indexing:

```
nextflow run singleron-RD/scrna \
 --star_genome {path to generated STAR genome index folder} \
 ...
```

> [!NOTE] > `--fasta` and `--gtf` are not required if `--star_genome` is provided.

## Input

samplesheet.csv:

```
sample,fastq_1,fastq_2
Sample_X,https://github.com/singleron-RD/scrna_test_data/raw/master/GEXSCOPE-V2/test_R1.fastq.gz,https://github.com/singleron-RD/scrna_test_data/raw/master/GEXSCOPE-V2/test_R2.fastq.gz
Sample_Y,https://github.com/singleron-RD/scrna_test_data/raw/master/GEXSCOPE-V2/Sample_Y_S1_L001_R1_001.fastq.gz,https://github.com/singleron-RD/scrna_test_data/raw/master/GEXSCOPE-V2/Sample_Y_S1_L001_R2_001.fastq.gz
Sample_Y,https://github.com/singleron-RD/scrna_test_data/raw/master/GEXSCOPE-V2/Sample_Y_S1_L002_R1_001.fastq.gz,https://github.com/singleron-RD/scrna_test_data/raw/master/GEXSCOPE-V2/Sample_Y_S1_L002_R2_001.fastq.gz
```

Each row represents a pair of fastq files.

## Main Output

Genome

- `genome/{genome_name}/` STAR genome index folder.

Matrix

- `matrix/{sample}.matrix/raw` Gene expression matrix file contains all barcodes(background + cell) from the barcode whitelist.

- `matrix/{sample}.matrix/filtered` Gene expression matrix file contains only cell barcodes. This file should be used as input to downstream analysis tools such as Seurat and Scanpy.

BAM

- `bam/{sample}.Aligned.sortedByCoord.out.bam` This bam file contains coordinate-sorted reads aligned to the genome.

Report

- `multiqc/multiqc_report.html` HTML QC report for reads and cells.
