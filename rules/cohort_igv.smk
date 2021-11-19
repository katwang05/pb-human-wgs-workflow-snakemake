if singleton:
    small_variants = f"samples/{samples[0]}/whatshap/{samples[0]}.{ref}.deepvariant.phased.vcf.gz"
    structural_variants = f"samples/{samples[0]}/pbsv/{samples[0]}.{ref}.pbsv.vcf.gz"
else:
    small_variants = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.glnexus.phased.vcf.gz"
    structural_variants = f"cohorts/{cohort}/pbsv/{cohort}.{ref}.pbsv.vcf.gz"


rule generate_igv_session:
    input:
        small_variants = small_variants,
        small_variants_tbi = [f"{x}.tbi" for x in small_variants],
        rare_small_variants = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf.gz",
        rare_small_variants_tbi = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf.gz.tbi"
        structural_variants = structural_variants,
        structural_variants_tbi = [f"{x}.tbi" for x in structural_variants],
        rare_structural_variants = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.vcf.gz",
        rare_structural_variants_tbi = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.vcf.gz.tbi"
        abams = [f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam" for sample in samples],
        abams_bai = [f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam.bai" for sample in samples],
        asm_bams = [f"samples/{sample}/hifiasm/{sample}.asm.{ref}.bam" for sample in samples],
        asm_bams_bai = [f"samples/{sample}/hifiasm/{sample}.asm.{ref}.bam.bai" for sample in samples],
        segdups = config['ref']['segdups'],
        oddregions = config['ref']['oddregions'],
        repeats = config['ref']['repeats']
    output:
        igv_session = f"cohorts/{cohort}/{cohort}.{ref}.igv.xml",
    log: f"cohorts/{cohort}/logs/igv/{cohort}.{ref}.igv.log"
    benchmark: f"cohorts/{cohort}/benchmarks/igv/{cohort}.{ref}.igv.tsv"
    conda: "envs/jinja.yaml"
    shell:
        f"""
        (python3 workflow/scripts/generate_igv_session.py \
            --webroot {config['webroot']} \
            --reference {ref} \
            --cohort {cohort} \
            --samples {','.join(samples)} \
            --segdups {{input.segdups}} \
            --oddregions {{input.oddregions}} \
            --repeats {{input.repeats}} \
            --output {{output.igv_session}}) 2>&1 {{log}}
        """
