rule cohort_pairs:
    input:
        vcf1 = "cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf.gz",
        vcf2 = lambda wildcards: f"cohorts/{cohort}/slivar/{pairs_dict[wildcards.sample]['pair']}.{ref}.deepvariant.phased.slivar.vcf.gz
    output:
        diff = "cohorts/{cohort}/diff/{cohort}.{ref}.diff.vcf"
    shell:
        "vcftools --vcf {input.vcf1} --diff {input.vcf2} --out {output.diff}"
