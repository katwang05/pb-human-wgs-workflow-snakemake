rule svpack_pairs:
    input:
        vcf1 = lambda wildcards: f"cohorts/{cohort}/svpack/{pairs_dict[wildcards.sample]['pair1']}.{ref}.pbsv.svpack.vcf.gz",
        vcf2 = lambda wildcards: f"cohorts/{cohort}/svpack/{pairs_dict[wildcards.sample]['pair2']}.{ref}.pbsv.svpack.vcf.gz"
    output:
        diff = "svpack/cohorts/{cohort}.pairs.{ref}.diff.vcf.gz"
    shell:
        "vcftools --vcf {input.vcf1} --diff {input.vcf2} --out {output.diff}"

rule slivar_pairs:
    input:
        vcf1 = lambda wildcards: f"cohorts/{cohort}/slivar/{pairs_dict[wildcards.sample]['pair1']}.{ref}.deepvariant.phased.slivar.vcf.gz",
        vcf2 = lambda wildcards: f"cohorts/{cohort}/slivar/{pairs_dict[wildcards.sample]['pair2']}.{ref}.deepvariant.phased.slivar.vcf.gz"
    output:
        diff = "slivar/cohorts/{cohort}.pairs.{ref}.diff.vcf.gz"
    shell:
        "vcftools --vcf {input.vcf1} --diff {input.vcf2} --out {output.diff}"
        
