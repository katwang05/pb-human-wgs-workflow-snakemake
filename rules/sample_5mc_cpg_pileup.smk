rule cpg_pileup:
    input:
        bam = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam",
        bai = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam.bai",
        reference = config['ref']['fasta'],
    output:
        [f"samples/{sample}/5mc_cpg_pileup/{sample}.{ref}.{infix}.{suffix}"
         for infix in ['combined', 'hap1', 'hap2']
         for suffix in ['bed', 'bw']]
    log: f"samples/{sample}/logs/5mc_cpg_pileup/{sample}.{ref}.log"
    benchmark: f"samples/{sample}/benchmarks/5mc_cpg_pileup/{sample}.{ref}.tsv"
    params:
        min_mapq = 1,
        min_coverage = 10,
        model = "pileup_calling_model.v1.tflite",
        prefix = f"samples/{sample}/5mc_cpg_pileup/{sample}.{ref}"
    threads: 8
    container: "docker://quay.io/pacbio/pb-cpg-tools:v2.1.0"
    shell:
        """
        (aligned_bam_to_cpg_scores \
			--threads {threads} \
			--bam {input.bam} \
			--ref {input.reference} \
			--output-prefix {params.prefix} \
			--min-mapq {params.min_mapq} \
			--min-coverage {params.min_coverage} \
			--model $PILEUP_MODEL_DIR/{params.model} > {log} 2>&1)
        """
