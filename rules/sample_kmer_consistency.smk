rule check_kmer_consistency:
    input:
        ref_modimers = config['ref']['modimers'],
        movie_modimers = expand(f"samples/{sample}/jellyfish/{{movie}}.modimers.tsv.gz", movie=movies)
    output: f"samples/{sample}/jellyfish/{sample}.kmerconsistency.txt"
    log: f"samples/{sample}/logs/jellyfish/kmerconsistency/{sample}.log"
    benchmark: f"samples/{sample}/benchmarks/jellyfish/kmerconsistency/{sample}.tsv"
    conda: "envs/base_python3.yaml"
    shell: "(python3 workflow/scripts/check_kmer_consistency.py {input.ref_modimers} {input.movie_modimers} > {output}) > {log} 2>&1"
