localrules: bgzip_vcf, tabix_vcf, tabix_bcf, create_ped, calculate_phrank


rule bgzip_vcf:
    input: f"cohorts/{cohort}/{{prefix}}.vcf"
    output: f"cohorts/{cohort}/{{prefix}}.vcf.gz"
    log: f"cohorts/{cohort}/logs/bgzip/{{prefix}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/bgzip/{{prefix}}.tsv"
    threads: 2
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Compressing {input}."
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule tabix_vcf:
    input: f"cohorts/{cohort}/{{prefix}}.vcf.gz"
    output: f"cohorts/{cohort}/{{prefix}}.vcf.gz.tbi"
    log: f"cohorts/{cohort}/logs/tabix/index/{{prefix}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/tabix/index/{{prefix}}.tsv"
    params: "-p vcf"
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Indexing {input}."
    shell: "(tabix {params} {input}) > {log} 2>&1"


rule tabix_bcf:
    input: f"cohorts/{cohort}/{{prefix}}.bcf"
    output: temp(f"cohorts/{cohort}/{{prefix}}.bcf.csi")
    log: f"cohorts/{cohort}/logs/tabix/index/{{prefix}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/tabix/index/{{prefix}}.tsv"
    params: "-p bcf"
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Indexing {input}."
    shell: "(tabix {params} {input}) > {log} 2>&1"


rule create_ped:
    input: config['cohort_yaml']
    output: f"cohorts/{cohort}/{cohort}.ped"
    log: f"cohorts/{cohort}/logs/yaml2ped/{cohort}.log"
    conda: "envs/pyyaml.yaml"
    message: f"Executing {{rule}}: Creating pedigree file for {cohort}."
    shell: f"(python3 workflow/scripts/yaml2ped.py {{input}} {cohort} {{output}}) > {{log}} 2>&1"


rule calculate_phrank:
    input:
        hpoterms = config['hpo']['terms'],
        hpodag = config['hpo']['dag'],
        hpoannotations = config['hpo']['annotations'],
        ensembltohgnc = config['ensembl_to_hgnc'],
        cohort_yaml = config['cohort_yaml']
    output: f"cohorts/{cohort}/{cohort}_phrank.tsv"
    log: f"cohorts/{cohort}/logs/calculate_phrank/{cohort}.log"
    conda: "envs/pyyaml.yaml"
    message: f"Executing {{rule}}: Calculate Phrank scores for {cohort}."
    shell:
        f"""(python3 workflow/scripts/calculate_phrank.py \
        {{input.hpoterms}} {{input.hpodag}} {{input.hpoannotations}} \
        {{input.ensembltohgnc}} {{input.cohort_yaml}} {cohort} {{output}}) > {{log}} 2>&1
        """


rule reformat_ensembl_gff:
    output: config['ref']['ensembl_gff']
    log: "logs/reformat_ensemble_gff.log"
    params: url = config['ref']['ensembl_gff_url']
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Downloaded and reformatting ensembl GFF to {output}."
    shell:
        """
        (wget -qO - {params.url} | zcat \
            | awk -v OFS="\t" '{{ if ($1=="##sequence-region") && ($2~/^G|K/) {{ print $0; }} else if ($0!~/G|K/) {{ print "chr" $0; }} }}' \
            | bgzip > {output}) > {log} 2>&1
        """


rule generate_lof_lookup:
    output: config['lof_lookup']
    log: "logs/generate_lof_lookup.log"
    params: url = config['lof_lookup_url']
    message: "Executing {rule}: Generating a lookup table for loss-of-function metrics at {output}."
    shell:
        """
        (wget -qO - {params.url} | zcat | cut -f 1,21,24 | tail -n+2 \
            | awk "{{ printf(\\"%s\\tpLI=%.3g;oe_lof=%.5g\\n\\", \$1, \$2, \$3) }}" > {output}) > {log} 2>&1
        """


rule generate_clinvar_lookup:
    output: config['clinvar_lookup']
    log: "logs/generate_clinvar_lookup.log"
    params: url = config['clinvar_lookup_url']
    message: "Executing {rule}: Generating a lookup table for clinvar gene descriptions at {output}."
    shell: "(wget -qO - {params.url} | cut -f 2,5 | grep -v ^$'\t' > {output}) > {log} 2>&1"
