# gvcf_list -> list of g.vcf.gz for all samples

rule glnexus:
    input:
        gvcf = gvcf_list,
        tbi = [f"{x}.tbi" for x in gvcf_list]
    output:
        bcf = temp(f"cohorts/{cohort}/glnexus/{cohort}.{ref}.deepvariant.glnexus.bcf"),
        scratch_dir = temp(directory(f"cohorts/{cohort}/glnexus/{cohort}.{ref}.GLnexus.DB/"))
    log: f"cohorts/{cohort}/logs/glnexus/{cohort}.{ref}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/glnexus/{cohort}.{ref}.tsv"
    container: f"docker://ghcr.io/dnanexus-rnd/glnexus:{config['GLNEXUS_VERSION']}"
    threads: 24
    message: f"Executing {{rule}}: Joint calling variants from {cohort} cohort."
    shell:
        """
        (rm -rf {output.scratch_dir} && \
        glnexus_cli --threads {threads} \
            --dir {output.scratch_dir} \
            --config DeepVariant_unfiltered {input.gvcf} > {output.bcf}) 2> {log}
        """


rule bcftools_bcf2vcf:
    input: f"cohorts/{cohort}/glnexus/{{prefix}}.bcf"
    output: f"cohorts/{cohort}/glnexus/{{prefix}}.vcf.gz"
    log: f"cohorts/{cohort}/logs/bcftools/view/glnexus/{{prefix}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/bcftools/view/glnexus/{{prefix}}.tsv"
    params: "--threads 4 -Oz"
    threads: 4
    conda: "envs/bcftools.yaml"
    message: "Executing {rule}: Converting GLnexus BCF to VCF for {input}."
    shell: "(bcftools view {params} {input} -o {output}) > {log} 2>&1"


rule split_glnexus_vcf:
    input:
        vcf = f"cohorts/{cohort}/glnexus/{cohort}.{ref}.deepvariant.glnexus.vcf.gz",
        tbi = f"cohorts/{cohort}/glnexus/{cohort}.{ref}.deepvariant.glnexus.vcf.gz.tbi"
    output: temp(f"cohorts/{cohort}/whatshap/regions/{cohort}.{ref}.{{region}}.deepvariant.glnexus.vcf")
    log: f"cohorts/{cohort}/logs/tabix/query/{cohort}.{ref}.{{region}}.glnexus.vcf.log"
    benchmark: f"cohorts/{cohort}/benchmarks/tabix/query/{cohort}.{ref}.{{region}}.glnexus.vcf.tsv"
    params: region = lambda wildcards: wildcards.region, extra = '-h'
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Extracting {wildcards.region} variants from {input}."
    shell: "(tabix {params.extra} {input.vcf} {params.region} > {output}) 2> {log}"


rule whatshap_phase:
    input:
        reference = config['ref']['fasta'],
        vcf = f"cohorts/{cohort}/whatshap/regions/{cohort}.{ref}.{{chromosome}}.deepvariant.glnexus.vcf.gz",
        tbi = f"cohorts/{cohort}/whatshap/regions/{cohort}.{ref}.{{chromosome}}.deepvariant.glnexus.vcf.gz.tbi",
        phaseinput = abam_list,
        phaseinputindex = [f"{x}.bai" for x in abam_list]
    output: temp(f"cohorts/{cohort}/whatshap/regions/{cohort}.{ref}.{{chromosome}}.deepvariant.glnexus.phased.vcf.gz")
    log: f"cohorts/{cohort}/logs/whatshap/phase/{cohort}.{ref}.{{chromosome}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/whatshap/phase/{cohort}.{ref}.{{chromosome}}.tsv"
    params:
        chromosome = lambda wildcards: wildcards.chromosome,
        extra = "--indels"
    conda: "envs/whatshap.yaml"
    message: "Executing {rule}: Phasing {input.vcf} using {input.phaseinput} for chromosome {wildcards.chromosome}."
    shell:
        """
        (whatshap phase {params.extra} \
            --chromosome {wildcards.chromosome} \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} {input.phaseinput}) > {log} 2>&1
        """

## to include pedigree information in phasing
# extra = """--indels --ped cohorts/{cohort}/{cohort}.ped \
#            --no-genetic-haplotyping --recombination-list cohorts/{cohort}/whatshap/{cohort}.recombination.list"""


rule whatshap_bcftools_concat:
    input:
        calls = expand(f"cohorts/{cohort}/whatshap/regions/{cohort}.{ref}.{{region}}.deepvariant.glnexus.phased.vcf.gz", region=all_chroms),
        indices = expand(f"cohorts/{cohort}/whatshap/regions/{cohort}.{ref}.{{region}}.deepvariant.glnexus.phased.vcf.gz.tbi", region=all_chroms)
    output: f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.glnexus.phased.vcf.gz"
    log: f"cohorts/{cohort}/logs/bcftools/concat/{cohort}.{ref}.whatshap.log"
    benchmark: f"cohorts/{cohort}/benchmarks/bcftools/concat/{cohort}.{ref}.whatshap.tsv"
    params: "-a -Oz"
    conda: "envs/bcftools.yaml"
    message: "Executing {rule}: Concatenating WhatsHap phased VCFs: {input.calls}"
    shell: "(bcftools concat {params} -o {output} {input.calls}) > {log} 2>&1"


rule whatshap_stats:
    input:
        vcf = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.glnexus.phased.vcf.gz",
        tbi = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.glnexus.phased.vcf.gz.tbi",
        chr_lengths = config['ref']['chr_lengths']
    output:
        gtf = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.glnexus.phased.gtf",
        tsv = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.glnexus.phased.tsv",
        blocklist = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.glnexus.phased.blocklist"
    log: f"cohorts/{cohort}/logs/whatshap/stats/{cohort}.{ref}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/whatshap/stats/{cohort}.{ref}.tsv"
    conda: "envs/whatshap.yaml"
    message: "Executing {rule}: Calculating phasing stats for {input.vcf}."
    shell:
        """
        (whatshap stats \
            --gtf {output.gtf} \
            --tsv {output.tsv} \
            --block-list {output.blocklist} \
            --chr-lengths {input.chr_lengths} \
            {input.vcf}) > {log} 2>&1
        """


# TODO: cleanup whatshap intermediates
