rule cpg_pileup:
    input:
        bam = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam",
        bai = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam.bai",
        reference = config['ref']['fasta'],
    output:
        [f"samples/{sample}/cpg_pileup/{sample}.{ref}.{infix}.denovo.{suffix}"
         for infix in ['combined', 'hap1', 'hap2']
         for suffix in ['bed', 'bw', 'mincov10.bed', 'mincov10.bw']]
    log: f"samples/{sample}/logs/cpg_pileup/{sample}.{ref}.log"
    benchmark: f"samples/{sample}/benchmarks/cpg_pileup/{sample}.{ref}.tsv"
    params:
        min_mapq = 1,
        min_coverage = 10,
        pileup_mode = "model",
        model_dir = "workflow/scripts/pb-CpG-tools/pileup_calling_model",
        modsites = "denovo",
        prefix = f"samples/{sample}/cpg_pileup/{sample}"
    threads: 48
    conda: "envs/pb-cpg-tools.yaml"
    shell:
        """
        (python3 workflow/scripts/pb-CpG-tools/aligned_bam_to_cpg_scores.py \
            -b {input.bam} -f {input.reference} -o {params.prefix} \
            -t {threads} -q {params.min_mapq} -m {params.modsites} \
            -p {params.pileup_mode} -d {params.model_dir} -c {params.min_coverage}) > {log} 2>&1
        """
