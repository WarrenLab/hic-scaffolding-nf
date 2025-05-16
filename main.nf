#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.extraYahsArgs = ''

process PRINT_VERSIONS {
    output:
    path("versions.txt")

    """
    echo "Chromap: \$(chromap --version 2>&1)" > versions.txt
    echo "YAHS: \$(yahs --version)" >> versions.txt
    java -jar $params.juicerToolsJar -V | grep Version >> versions.txt
    echo "assembly-stats: \$(assembly-stats -v)" >> versions.txt
    """
}

process SAMTOOLS_FAIDX {
    input:
    path(contigsFasta)

    output:
    path("${contigsFasta}.fai")

    """
    samtools faidx $contigsFasta
    """
}

process CHROMAP_INDEX {
    input:
    path(contigsFasta)

    output:
    path("contigs.index")

    """
    chromap -i -r $contigsFasta -o contigs.index
    """
}

process CHROMAP_ALIGN {
    input:
    path(contigsFasta)
    path(contigsChromapIndex)
    path(r1Reads)
    path(r2Reads)

    output:
    path("aligned.bam")

    """
    chromap \
        --preset hic \
        -r $contigsFasta \
        -x $contigsChromapIndex \
        --remove-pcr-duplicates \
        -1 $r1Reads \
        -2 $r2Reads \
        --SAM \
        -o aligned.sam \
        -t ${task.cpus}

    samtools view -bh aligned.sam | samtools sort -n > aligned.bam
    """
}

process YAHS_SCAFFOLD {
    input:
    path("contigs.fa")
    path("contigs.fa.fai")
    path("aligned.bam")

    output:
    path("yahs.out.bin"), emit: bin
    path("yahs.out_scaffolds_final.agp"), emit: agp
    path("yahs.out_scaffolds_final.fa"), emit: fasta

    """
    yahs $params.extraYahsArgs contigs.fa aligned.bam
    """
}

#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.extraYahsArgs = ''

process PRINT_VERSIONS {
    output:
    path("versions.txt")

    """
    echo "Chromap: \$(chromap --version 2>&1)" > versions.txt
    echo "YAHS: \$(yahs --version)" >> versions.txt
    java -jar $params.juicerToolsJar -V | grep Version >> versions.txt
    echo "assembly-stats: \$(assembly-stats -v)" >> versions.txt
    """
}

process SAMTOOLS_FAIDX {
    input:
    path(contigsFasta)

    output:
    path("${contigsFasta}.fai")

    """
    samtools faidx $contigsFasta
    """
}

process CHROMAP_INDEX {
    input:
    path(contigsFasta)

    output:
    path("contigs.index")

    """
    chromap -i -r $contigsFasta -o contigs.index
    """
}

process CHROMAP_ALIGN {
    input:
    path(contigsFasta)
    path(contigsChromapIndex)
    path(r1Reads)
    path(r2Reads)

    output:
    path("aligned.bam")

    """
    chromap \
        --preset hic \
        -r $contigsFasta \
        -x $contigsChromapIndex \
        --remove-pcr-duplicates \
        -1 $r1Reads \
        -2 $r2Reads \
        --SAM \
        -o aligned.sam \
        -t ${task.cpus}

    samtools view -bh aligned.sam | samtools sort -n > aligned.bam
    """
}

process YAHS_SCAFFOLD {
    input:
    path("contigs.fa")
    path("contigs.fa.fai")
    path("aligned.bam")

    output:
    path("yahs.out.bin"), emit: bin
    path("yahs.out_scaffolds_final.agp"), emit: agp
    path("yahs.out_scaffolds_final.fa"), emit: fasta

    """
    yahs $params.extraYahsArgs contigs.fa aligned.bam
    """
}

process JUICER_PRE {
    input:
    path("yahs.out.bin")
    path("yahs.out_scaffolds_final.agp")
    path("contigs.fa.fai")

    output:
    path("manual_JBAT.*")
    path("manual_JBAT.hic")

    script:
    """
    # Run juicer pre, redirect output to log
    juicer pre -a -o manual_JBAT yahs.out.bin yahs.out_scaffolds_final.agp contigs.fa.fai > manual_JBAT.log 2>&1

    # Extract chromosome sizes from the log file and pass to juicer_tools pre
    java -Xmx32G -jar $params.juicerToolsJar pre manual_JBAT.txt manual_JBAT.hic <(grep PRE_C_SIZE manual_JBAT.log | awk '{print \$2, \$3}') && \
    mv manual_JBAT.hic.part manual_JBAT.hic
    """
}
process ASSEMBLY_STATS {
    input:
    path("yahs.out_scaffolds_final.fa")

    output:
    path("assembly_stats.txt")

    """
    assembly-stats yahs.out_scaffolds_final.fa > assembly_stats.txt
    """
}

workflow {
    // TODO do a parameter check
    PRINT_VERSIONS()

    r1Reads = Channel.fromPath(params.r1Reads)
    r2Reads = Channel.fromPath(params.r2Reads)
    contigs = Channel.fromPath(params.contigs)

    SAMTOOLS_FAIDX(contigs)
    CHROMAP_INDEX(contigs)

    CHROMAP_ALIGN(contigs, CHROMAP_INDEX.out, r1Reads, r2Reads)

    YAHS_SCAFFOLD(contigs, SAMTOOLS_FAIDX.out, CHROMAP_ALIGN.out)

    JUICER_PRE(YAHS_SCAFFOLD.out.bin, YAHS_SCAFFOLD.out.agp, SAMTOOLS_FAIDX.out)

    ASSEMBLY_STATS(YAHS_SCAFFOLD.out.fasta)
}

process ASSEMBLY_STATS {
    input:
    path("yahs.out_scaffolds_final.fa")

    output:
    path("assembly_stats.txt")

    """
    assembly-stats yahs.out_scaffolds_final.fa > assembly_stats.txt
    """
}

workflow {
    // TODO do a parameter check
    PRINT_VERSIONS()

    r1Reads = Channel.fromPath(params.r1Reads)
    r2Reads = Channel.fromPath(params.r2Reads)
    contigs = Channel.fromPath(params.contigs)

    SAMTOOLS_FAIDX(contigs)
    CHROMAP_INDEX(contigs)

    CHROMAP_ALIGN(contigs, CHROMAP_INDEX.out, r1Reads, r2Reads)

    YAHS_SCAFFOLD(contigs, SAMTOOLS_FAIDX.out, CHROMAP_ALIGN.out)

    JUICER_PRE(YAHS_SCAFFOLD.out.bin, YAHS_SCAFFOLD.out.agp, SAMTOOLS_FAIDX.out)

    ASSEMBLY_STATS(YAHS_SCAFFOLD.out.fasta)
}

