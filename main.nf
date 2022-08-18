#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process SAMTOOLS_FAIDX {
    input:
    file(contigsFasta)

    output:
    file("${contigsFasta}.fai")

    """
    samtools index $contigsFasta
    """
}

process CHROMAP_INDEX {
    input:
    file(contigsFasta)

    output:
    file("contigs.index")

    """
    chromap -i -r $contigsFasta -o contigs.index
    """
}

process CHROMAP_ALIGN {
    input:
    file(referenceFasta)
    file(referenceChromapIndex)
    file(r1Reads)
    file(r2Reads)

    output:
    file("aligned.bam")

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
    file("contigs.fa")
    file("contigs.fa.fai")
    file("aligned.bam")

    output:
    file("yahs.out.bin"), emit: bin
    file("yahs.out_scaffolds_final.agp"), emit: agp
    file("yahs.out_scaffolds_final.fa"), emit: fasta

    """
    yahs contigs.fa aligned.bam
    """
}

process JUICER_PRE {
    input:
    file("yahs.out.bin")
    file("yahs.out_scaffolds_final.agp")
    file("contigs.fa.fai")

    output:
    file("out_JBAT.*")

    // TODO get asm size
    """
    juicer pre -a -o out_JBAT \
        yahs.out.bin \
        yahs.out_scaffolds_final.agp \
        contigs.fa.fai

    java -Xmx36G -jar $params.juicerToolsJar \
        pre out_JBAT.txt out_JBAT.hic <(echo "assembly \${asm_size}")
    """
}

workflow {
    // TODO do a parameter check

    r1Reads = Channel.fromPath(params.r1Reads)
    r2Reads = Channel.fromPath(params.r2Reads)
    contigs = Channel.fromPath(params.contigs)

    SAMTOOLS_FAIDX(contigs)
    CHROMAP_INDEX(contigs)

    CHROMAP_ALIGN(contigs, CHROMAP_INDEX.out, r1Reads, r2Reads)

    YAHS_SCAFFOLD(contigs, SAMTOOLS_FAIDX.out, CHROMAP_ALIGN.out)

    JUICER_PRE(YAHS_SCAFFOLD.out.bin, YAHS_SCAFFOLD.out.agp, SAMTOOLS_FAIDX.out)
}
