rule jellyfish_merge:
    input: expand(f"samples/{sample}/jellyfish/{{movie}}.jf", movie=movies)
    output: f"samples/{sample}/jellyfish/{sample}.jf"
    log: f"samples/{sample}/logs/jellyfish/merge/{sample}.log"
    benchmark: f"samples/{sample}/benchmarks/jellyfish/merge/{sample}.tsv"
    conda: "envs/jellyfish.yaml"
    shell:
        f"""
        if [[ "{{input}}" =~ " " ]]
        then
            (jellyfish merge -o {{output}} {{input}}) > {{log}} 2>&1
        else
            ln -s {movies[0]}.jf {{output}}
        fi
        """
