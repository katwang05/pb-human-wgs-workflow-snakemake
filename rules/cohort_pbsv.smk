# svsig_dict[chrom] = list of svsigs for all movies for all samples for given chrom

localrules: bcftools_concat_pbsv_vcf


rule pbsv_call:
    input:
        svsigs = lambda wildcards: svsig_dict[wildcards.region],
        reference = config['ref']['fasta']
    output: temp(f"cohorts/{cohort}/pbsv/{cohort}.{ref}.chrom_vcfs/{cohort}.{ref}.{{region}}.pbsv.vcf")
    log: f"cohorts/{cohort}/logs/pbsv/call/{cohort}.{ref}.{{region}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/pbsv/call/{cohort}.{ref}.{{region}}.tsv"
    params:
        region = lambda wildcards: wildcards.region,
        extra = "--hifi -m 20",
        loglevel = "INFO"
    threads: 8
    conda: "envs/pbsv.yaml"
    message: "Executing {rule}: Calling structural variants in {wildcards.region} from SVSIGs: {input.svsigs}"
    shell:
        """
        (pbsv call {params.extra} \
            --log-level {params.loglevel} \
            --num-threads {threads} \
            {input.reference} {input.svsigs} {output}) > {log} 2>&1
        """


rule bcftools_concat_pbsv_vcf:
    input:
        calls = expand(f"cohorts/{cohort}/pbsv/{cohort}.{ref}.chrom_vcfs/{cohort}.{ref}.{{region}}.pbsv.vcf.gz", region=all_chroms),
        indices = expand(f"cohorts/{cohort}/pbsv/{cohort}.{ref}.chrom_vcfs/{cohort}.{ref}.{{region}}.pbsv.vcf.gz.tbi", region=all_chroms)
    output: f"cohorts/{cohort}/pbsv/{cohort}.{ref}.pbsv.vcf"
    log: f"cohorts/{cohort}/logs/bcftools/concat/{cohort}.{ref}.pbsv.vcf.log"
    benchmark: f"cohorts/{cohort}/benchmarks/bcftools/concat/{cohort}.{ref}.pbsv.vcf.tsv"
    conda: "envs/bcftools.yaml"
    message: "Executing {rule}: Concatenating pbsv VCFs: {input.calls}"
    shell: "(bcftools concat -a -o {output} {input.calls}) > {log} 2>&1"
