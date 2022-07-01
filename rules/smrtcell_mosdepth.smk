rule mosdepth:
    input:
        bam = f"samples/{{sample}}/aligned/{{prefix}}.{ref}.bam",
        bai = f"samples/{{sample}}/aligned/{{prefix}}.{ref}.bam.bai"
    output:
        f"samples/{{sample}}/mosdepth/{{prefix}}.{ref}.mosdepth.global.dist.txt",
        f"samples/{{sample}}/mosdepth/{{prefix}}.{ref}.mosdepth.region.dist.txt",
        f"samples/{{sample}}/mosdepth/{{prefix}}.{ref}.mosdepth.summary.txt",
        f"samples/{{sample}}/mosdepth/{{prefix}}.{ref}.regions.bed.gz"
    log: f"samples/{{sample}}/logs/mosdepth/{{prefix}}.{ref}.log"
    benchmark: f"samples/{{sample}}/benchmarks/mosdepth/{{prefix}}.{ref}.tsv"
    params:
        by = "500",
        prefix = f"samples/{{sample}}/mosdepth/{{prefix}}.{ref}",
        extra = "--no-per-base --use-median"
    threads: 4
    conda: "envs/mosdepth.yaml"
    shell:
        """
        (mosdepth \
            --threads {threads} --by {params.by} \
            {params.extra} {params.prefix} {input.bam}) > {log} 2>&1
        """
