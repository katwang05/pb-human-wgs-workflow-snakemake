rule bcftools_concat_compound_het_vcfs:
    input:
        small_vcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf.gz",
        small_tbi = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf.gz.tbi",
        sv_vcf = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.slivar.vcf.gz",
        sv_tbi = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.slivar.vcf.gz.tbi",
    output:
        temp(f"cohorts/{cohort}/slivar/{cohort}.{ref}.concat.vcf"),
    log:
        f"cohorts/{cohort}/bcftools/concat/{cohort}.{ref}.log"
    benchmark:
        f"cohorts/{cohort}/bcftools/concat/{cohort}.{ref}.tsv"
    params:
    conda:
        "envs/bcftools.yaml"
    message: f"Executing {{rule}}: Merging small variant and structural variant VCFs for cohort {cohort}."
    shell: "(bcftools concat -a -o {output} {input.small_vcf} {input.sv_vcf}) > {log} 2>&1"


skip_list = [
    'non_coding_transcript',
    'intron',
    'non_coding_transcript',
    'non_coding',
    'upstream_gene',
    'downstream_gene',
    'non_coding_transcript_exon',
    'NMD_transcript',
    '5_prime_UTR',
    '3_prime_UTR',
    'sv:intron',
    'sv:utr'
    ]


rule all_compound_hets:
    input: 
        vcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.concat.vcf.gz",
        tbi = f"cohorts/{cohort}/slivar/{cohort}.{ref}.concat.vcf.gz.tbi",
        ped = f"cohorts/{cohort}/{cohort}.ped"
    output:
        vcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.concat.compound-hets.vcf"
    log: f"cohorts/{cohort}/logs/slivar/compound-hets/{cohort}.{ref}.concat.compound-hets.log"
    benchmark: f"cohorts/{cohort}/benchmarks/slivar/compound-hets/{cohort}.{ref}.concat.compound-hets.tsv"
    params: skip = ",".join(skip_list)
    conda: "envs/slivar.yaml"
    message: f"Executing {{rule}}: Finding compound hets in {cohort}."
    shell:
        """
        (slivar compound-hets \
            --skip {params.skip} \
            --vcf {input.vcf} \
            --sample-field comphet_side \
            --ped {input.ped} \
            --allow-non-trios \
            > {output.vcf}) > {log} 2>&1
        """


info_fields = [
    'gnomad_af',
    'hprc_af',
    'gnomad_nhomalt',
    'hprc_nhomalt',
    'gnomad_ac',
    'hprc_ac',
    'SVTYPE',
    'SVLEN',
    'SVANN',
    'CIPOS',
    'MATEID',
    'END'
    ]


rule all_compound_het_tsv:
    input:
        comphet_vcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.concat.compound-hets.vcf.gz",
        ped = f"cohorts/{cohort}/{cohort}.ped",
        lof_lookup = config['lof_lookup'],
        clinvar_lookup = config['clinvar_lookup'],
        phrank_lookup = f"cohorts/{cohort}/{cohort}_phrank.tsv"
    output:
        comphet_tsv = f"cohorts/{cohort}/slivar/{cohort}.{ref}.concat.compound-hets.tsv"
    log: f"cohorts/{cohort}/logs/slivar/tsv/{cohort}.{ref}.concat.compound-hets.log"
    benchmark: f"cohorts/{cohort}/benchmarks/slivar/tsv/{cohort}.{ref}.concat.compound-hets.tsv"
    params: info = "".join([f"--info-field {x} " for x in info_fields])
    conda: "envs/slivar.yaml"
    message: "Executing {rule}: Converting annotated VCFs to TSVs for easier interpretation."
    shell:
        """
        (slivar tsv \
            {params.info} \
            --sample-field slivar_comphet \
            --info-field slivar_comphet \
            --csq-field BCSQ \
            --gene-description {input.lof_lookup} \
            --gene-description {input.clinvar_lookup} \
            --gene-description {input.phrank_lookup} \
            --ped {input.ped} \
            --out /dev/stdout \
            {input.comphet_vcf} \
            | sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
            > {output.comphet_tsv}) > {log} 2>&1
        """
