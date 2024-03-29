## Introduction

**singleron-RD/scrna** is a bioinformatics pipeline for processing single-cell RNA-seq data. It is designed for analysis of data generated by the Singleron GEXSCOPE kit, but can be adapted to analyze single-cell RNAseq data from other platforms.
It takes a samplesheet and FASTQ files as input, performs the following steps:

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Auto-detect GEXSCOPE protocols.
3. Mapping, demultiplexing and quantification ([`STARSolo`](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md))
4. HTML report ([`MultiQC`](http://multiqc.info/))

## Documents

- [Quick Start](./docs/quickstart.md)
- [Usage](./docs/usage.md)
- [Parameters](./docs/parameters.md)
- [Output](./docs/output.md)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
