# hic-scaffolding-nf
Nextflow pipeline for scaffolding genome assemblies with Hi-C reads

## Introduction
This pipeline requires the following inputs:
1. A fasta file containing assembled contigs (`--contigs`)
2. Hi-C reads in paired-end fastq(.gz) format (`--r1Reads` and `--r2Reads`)

It then performs the following tasks:
1. Aligns the Hi-C reads to the contigs using [chromap][chromap]
2. Scaffolds the contigs using [yahs][yahs]
3. Prepares all the files you need to do manual curation in
   [Juicebox][juicer_tools]

and produces the following outputs:
* Alignments in bam format (`out/chromap/aligned.bam`)
* A scaffolded assembly in both agp and fasta formats
  (`out/scaffolds/yahs.out_scaffolds_final.[agp,fa]`)
* `.hic` and `.assembly` files for loading in Juicebox Assembly Tools
  (`out/juicebox_input/out_JBAT.[hic,assembly]`)

## Configuration
### Running on Lewis
If you're running this on the Lewis cluster, I've already got a profile set up
with everything you need, so just add `-profile lewis` to the command and
you're good to go.

### Running on another cluster/cloud/locally
This pipeline has the following dependencies:
* [nextflow][nextflow]
* [chromap][chromap]
* [yahs][yahs]
* [JuicerTools][juicer_tools]

Nextflow must be in your path. You can get nextflow to make a conda environment
containing chromap and yahs for you with `-profile conda` (note one dash!).
JuicerTools is distributed as a jar file, so you need to tell the pipeline
where it is by adding the argument `--juicer-tools-jar /path/to/jar` (note two
dashes!). You can also add this stuff to a config file called `nextflow.config`
in the directory from which you're running it (see nextflow documentation).

## Running
```bash
nextflow run WarrenLab/hic-scaffolding-nf \
    --contigs contigs.fa \
    --r1Reads hic_reads_R1.fastq.gz \
    --r2Reads hic_reads_R2.fastq.gz
```
You'll need to add a couple options depending on your configuration (see
section above).

[nextflow]: https://www.nextflow.io/
[chromap]: https://github.com/haowenz/chromap
[yahs]: https://github.com/c-zhou/yahs
[juicer_tools]: https://github.com/aidenlab/juicer/wiki/Download
