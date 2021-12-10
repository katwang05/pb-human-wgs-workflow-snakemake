rule svpack_filter_annotated:
    input:
        pbsv_vcf = svpack_input,
        eee_vcf = config['ref']['eee_vcf'],
        gnomadsv_vcf = config['ref']['gnomadsv_vcf'],
        hprc_pbsv_vcf = config['ref']['hprc_pbsv_vcf'],
        decode_vcf = config['ref']['decode_vcf'],
        gff = config['ref']['ensembl_gff']
    output: f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.vcf"
    log: f"cohorts/{cohort}/logs/svpack/{cohort}.{ref}.pbsv.svpack.log"
    benchmark: f"cohorts/{cohort}/benchmarks/svpack/{cohort}.{ref}.pbsv.svpack.tsv"
    params:
        min_sv_length = 50
    conda: "envs/svpack.yaml"
    shell:
        """
        (python workflow/scripts/svpack/svpack filter --pass-only \
            --min-svlen {params.min_sv_length} {input.pbsv_vcf} | \
            python workflow/scripts/svpack/svpack match -v - {input.eee_vcf} | \
            python workflow/scripts/svpack/svpack match -v - {input.gnomadsv_vcf} | \
            python workflow/scripts/svpack/svpack match -v - {input.hprc_pbsv_vcf} | \
            python workflow/scripts/svpack/svpack match -v - {input.decode_vcf} | \
            python workflow/scripts/svpack/svpack consequence - {input.gff} | \
            > {output}) 2> {log}
        """


if singleton:
    # singleton
    slivar_filters = [
        "--info 'variant.FILTER==\"PASS\"}'",
        "--family-expr 'recessive:fam.every(segregating_recessive)'",
        "--family-expr 'x_recessive:(variant.CHROM == \"chrX\") && fam.every(segregating_recessive_x)'",
        "--family-expr 'dominant:fam.every(segregating_dominant)'",
        "--family-expr 'x_dominant:(variant.CHROM == \"chrX\") && fam.every(segregating_dominant_x)'",
        "--sample-expr 'comphet_side:sample.het'"
        ]
else:
    # trio cohort
    slivar_filters = [
        "--info 'variant.FILTER==\"PASS\"'",
        "--family-expr 'recessive:fam.every(segregating_recessive)'",
        "--family-expr 'x_recessive:(variant.CHROM == \"chrX\") && fam.every(segregating_recessive_x)'",
        "--family-expr 'dominant:fam.every(segregating_dominant)'",
        "--family-expr 'x_dominant:(variant.CHROM == \"chrX\") && fam.every(segregating_dominant_x)'",
        "--trio 'comphet_side:comphet_side(kid, mom, dad) && kid.affected'"
    ]


rule slivar_structural_variant:
    input:
        vcf = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.vcf.gz",
        tbi = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.vcf.gz.tbi",
        ped = f"cohorts/{cohort}/{cohort}.ped",
        js = "workflow/scripts/slivar-functions-pbsv.js",
        ref = config['ref']['fasta']
    output: f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.slivar.vcf"
    log: f"cohorts/{cohort}/logs/slivar/filter/{cohort}.{ref}.pbsv.svpack.slivar.vcf.log"
    benchmark: f"cohorts/{cohort}/benchmarks/slivar/filter/{cohort}.{ref}.pbsv.svpack.slivar.tsv"
    params: filters = slivar_filters
    threads: 12
    conda: "envs/slivar.yaml"
    message: "Executing {rule}: Annotating {input.vcf} and applying filters."
    shell:
        """
        (pslivar --processes {threads} \
            --fasta {input.ref}\
            --pass-only \
            --js {input.js} \
            {params.filters} \
            --vcf {input.vcf} \
            --ped {input.ped} \
            --out-vcf {output}) > {log} 2>&1
        """


rule slivar_structural_variant_compound_hets:
    input: 
        vcf = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.slivar.vcf.gz",
        tbi = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.slivar.vcf.gz.tbi",
        ped = f"cohorts/{cohort}/{cohort}.ped"
    output:
        vcf = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.slivar.compound-hets.vcf"
    log: f"cohorts/{cohort}/logs/slivar/compound-hets/{cohort}.{ref}.pbsv.svpack.slivar.compound-hets.vcf.log"
    benchmark: f"cohorts/{cohort}/benchmarks/slivar/compound-hets/{cohort}.{ref}.pbsv.svpack.slivar.compound-hets.vcf.tsv"
    params: skip = ",".join(skip_list)
    conda: "envs/slivar.yaml"
    message: f"Executing {{rule}}: Finding compound hets in {cohort}."
    shell:
        """
        (slivar compound-hets \
            --vcf {input.vcf} \
            --sample-field comphet_side \
            --ped {input.ped} \
            --allow-non-trios \
            | python3 workflow/scripts/add_comphet_phase.py \
            > {output.vcf}) > {log} 2>&1
        """


info_fields = [
    'SVTYPE',
    'SVLEN',
    'SVANN',
    'CIPOS',
    'MATEID',
    'END'
    ]


rule slivar_svpack_tsv:
    input:
        filt_vcf = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.slivar.vcf.gz",
        comphet_vcf = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.slivar.compound-hets.vcf.gz"
        ped = f"cohorts/{cohort}/{cohort}.ped",
        lof_lookup = config['lof_lookup'],
        clinvar_lookup = config['clinvar_lookup'],
        phrank_lookup = f"cohorts/{cohort}/{cohort}_phrank.tsv"
    output:
        filt_tsv = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.slivar.tsv",
        comphet_tsv = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.slivar.compound-hets.tsv"
    log: f"cohorts/{cohort}/logs/slivar/tsv/{cohort}.{ref}.pbsv.svpack.log"
    benchmark: f"cohorts/{cohort}/benchmarks/slivar/tsv/{cohort}.{ref}.pbsv.svpack.tsv"
    params: info = "".join([f"--info-field {x} " for x in info_fields])
    conda: "envs/slivar.yaml"
    message: "Executing {rule}: Converting annotated VCFs to TSVs for easier interpretation."
    shell:
        """
        (slivar tsv \
            {params.info} \
            --sample-field dominant \
            --sample-field x_dominant \
            --sample-field recessive \
            --sample-field x_recessive \
            --csq-field BCSQ \
            --gene-description {input.lof_lookup} \
            --gene-description {input.clinvar_lookup} \
            --gene-description {input.phrank_lookup} \
            --ped {input.ped} \
            --out /dev/stdout \
            {input.filt_vcf} \
            | sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
            > {output.filt_tsv}
        slivar tsv \
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
