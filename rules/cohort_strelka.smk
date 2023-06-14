rule run_strelka2:
    input:
        tumor = "{sample}.tumor.bam",
        normal = "{sample}.normal.bam",
        ref = "reference.fasta",
        indels = "{sample}.indels.vcf"
    output:
        vcf = "{sample}.strelka.vcf"
    params:
        dir = "{sample}.strelka"
    shell:
        """
        strelka2 --tumor {input.tumor} --normal {input.normal} \
        --referenceFasta {input.ref} --indelCandidates {input.indels} \
        --runDir {params.dir}
        {params.dir}/runWorkflow.py
        """
