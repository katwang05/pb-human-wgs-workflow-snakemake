rule mosdepth:
    input:
        bam = f"samples/{{sample}}/whatshap/{{prefix}}.bam",
        bai = f"samples/{{sample}}/whatshap/{{prefix}}.bam.bai"
    output:
        f"samples/{{sample}}/mosdepth/{{prefix}}.mosdepth.global.dist.txt",
        f"samples/{{sample}}/mosdepth/{{prefix}}.mosdepth.region.dist.txt",
        f"samples/{{sample}}/mosdepth/{{prefix}}.mosdepth.summary.txt",
        f"samples/{{sample}}/mosdepth/{{prefix}}.regions.bed.gz"
    log: f"samples/{{sample}}/logs/mosdepth/{{prefix}}.log"
    benchmark: f"samples/{{sample}}/benchmarks/mosdepth/{{prefix}}.tsv"
    params:
        by = "500",
        prefix = f"samples/{{sample}}/mosdepth/{{prefix}}",
        extra = "--no-per-base --use-median"
    threads: 4
    conda: "envs/mosdepth.yaml"
    shell:
        """
        (mosdepth --threads {threads} --by {params.by} \
            {params.extra} {params.prefix} {input.bam}) > {log} 2>&1
        """
