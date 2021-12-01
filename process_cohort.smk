import os
import re
import yaml
from pathlib import Path
from collections import defaultdict


configfile: "workflow/reference.yaml"         # reference information
configfile: "workflow/config.yaml"            # general configuration
shell.prefix(f"set -o pipefail; umask 002; export TMPDIR={config['tmpdir']}; export SINGULARITY_TMPDIR={config['tmpdir']}; ")  # set g+w


def get_samples(cohortyaml=config['cohort_yaml'], cohort_id=config['cohort']):
    """Find all samples associated with cohort."""
    with open(cohortyaml, 'r') as yamlfile:
        cohort_list = yaml.load(yamlfile, Loader = yaml.FullLoader)
    for c in cohort_list:
        if c['id'] == cohort_id:
            break
    samples = []
    for affectedstatus in ['affecteds', 'unaffecteds']:
        if affectedstatus in c:
            for individual in range(len(c[affectedstatus])):
                samples.append(c[affectedstatus][individual]['id'])
    return samples

def get_trios(cohortyaml=config['cohort_yaml'], cohort_id=config['cohort']):
    """Find all trios associated with cohort."""
    with open(cohortyaml, 'r') as yamlfile:
        cohort_list = yaml.load(yamlfile, Loader = yaml.FullLoader)
    for c in cohort_list:
        if c['id'] == cohort_id:
            break
    trio_dict = defaultdict(dict)
    for affectedstatus in ['affecteds', 'unaffecteds']:
        if affectedstatus in c:
            for individual in range(len(c[affectedstatus])):
                if ('parents' in c[affectedstatus][individual]) \
                        and (len(c[affectedstatus][individual]['parents']) == 2):
                    trio_dict[c[affectedstatus][individual]['id']]['parent1'] = c[affectedstatus][individual]['parents'][0]
                    trio_dict[c[affectedstatus][individual]['id']]['parent2'] = c[affectedstatus][individual]['parents'][1]
    return trio_dict

# cohort will be provided at command line with `--config cohort=$COHORT`
cohort = config['cohort']
ref = config['ref']['shortname']
all_chroms = config['ref']['autosomes'] + config['ref']['sex_chrom'] + config['ref']['mit_chrom']
print(f"Processing cohort {cohort} with reference {ref}.")

# find all samples in cohort
samples = get_samples()

if len(samples) == 0:
    print(f"No samples in {cohort}.") and exit
elif len(samples) == 1:
    singleton = True
else:
    singleton = False
print(f"Samples in cohort: {samples}.")

# find all trios in cohort
trio_dict = get_trios()
if trio_dict:
    print(f"Parents listed for samples: {list(trio_dict.keys())}")
    for s in trio_dict.keys():
        print(f"Parent IDs for {s}: {trio_dict[s]['parent1']}, {trio_dict[s]['parent2']}")
else:
    print(f"No trios found in cohort.")

# scan smrtcells/ready directory for uBAMs or FASTQs that are ready to process
# uBAMs have priority over FASTQs in downstream processes if both are available
ubam_pattern = re.compile(r'smrtcells/ready/(?P<sample>[A-Za-z0-9_-]+)/(?P<movie>m\d{5}[Ue]?_\d{6}_\d{6}).(ccs|hifi_reads).bam')
ubam_dict = defaultdict(dict)
fastq_pattern = re.compile(r'smrtcells/ready/(?P<sample>[A-Za-z0-9_-]+)/(?P<movie>m\d{5}[Ue]?_\d{6}_\d{6}).fastq.gz')
fastq_dict = defaultdict(dict)
for infile in Path('smrtcells/ready').glob('**/*.bam'):
    ubam_match = ubam_pattern.search(str(infile))
    if ubam_match:
        # create a dict-of-dict to link samples to movie context to uBAM filenames
        ubam_dict[ubam_match.group('sample')][ubam_match.group('movie')] = str(infile)
for infile in Path('smrtcells/ready').glob('**/*.fastq.gz'):
    fastq_match = fastq_pattern.search(str(infile))
    if fastq_match:
        # create a dict-of-dict to link samples to movie context to FASTQ filenames
        fastq_dict[fastq_match.group('sample')][fastq_match.group('movie')] = str(infile)
# get list of unique movie names for each sample, ignoring redundancy between ubam and fastq
ubam_fastq_dict = {sample:list(set(list(ubam_dict[sample].keys()) + list(fastq_dict[sample].keys()))) for sample in list(ubam_dict.keys()) + list(fastq_dict.keys())}

# scan samples/*/aligned to generate a dict-of-lists-of-movies for 
pattern = re.compile(r'samples/(?P<sample>[A-Za-z0-9_-]+)/aligned/(?P<movie>m\d{5}[Ue]?_\d{6}_\d{6})\.(?P<reference>.*).bam')
movie_dict = defaultdict(list)
abam_list = []
for infile in Path(f"samples").glob('**/aligned/*.bam'):
    match = pattern.search(str(infile))
    if match and (match.group('sample') in samples) and (match.group('reference') == ref):
        movie_dict[match.group('sample')].append(match.group('movie'))
        abam_list.append(infile)

# singletons and cohorts provide different input to slivar and svpack
if singleton:
    # use the sample level VCFs
    slivar_input = f"samples/{samples[0]}/whatshap/{samples[0]}.{ref}.deepvariant.phased.vcf.gz"
    svpack_input = f"samples/{samples[0]}/pbsv/{samples[0]}.{ref}.pbsv.vcf.gz"
    gvcf_list = []   # unused
    svsig_dict = []  # unused
else:
    # generate joint-called VCFs
    slivar_input = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.glnexus.phased.vcf.gz"
    svpack_input = f"cohorts/{cohort}/pbsv/{cohort}.{ref}.pbsv.vcf.gz"
    gvcf_list = [f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.g.vcf.gz" for sample in samples]
    svsig_dict = {region: [f"samples/{sample}/pbsv/svsig/{movie}.{ref}.{region}.svsig.gz"
                           for sample in samples
                           for movie in movie_dict[sample]]
                  for region in all_chroms}


# build a list of targets
targets = []
include: 'rules/cohort_common.smk'

# assemble with hifiasm
include: 'rules/cohort_hifiasm.smk'
if 'trio_assembly' in config['cohort_targets']:
    # assembly and stats
    targets.extend([f"cohorts/{cohort}/hifiasm/{trio}.asm.dip.{infix}.{suffix}"
                for suffix in ['fasta.gz', 'fasta.stats.txt']
                for infix in ['hap1.p_ctg', 'hap2.p_ctg']
                for trio in trio_dict.keys()])
    # assembly alignments
    targets.extend([f"cohorts/{cohort}/hifiasm/{trio}.asm.{ref}.{suffix}"
                for suffix in ['bam', 'bam.bai']
                for trio in trio_dict.keys()])

# generate a cohort level pbsv vcf or use singleton vcf
include: 'rules/cohort_pbsv.smk'
if 'pbsv_vcf' in config['cohort_targets']:
    targets.extend([svpack_input, svpack_input + '.tbi'])

# annotate and filter pbsv vcf with svpack
include: 'rules/cohort_svpack.smk'
if 'svpack' in config['cohort_targets']:
    targets.extend([f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.{suffix}"
                    for suffix in ['vcf.gz', 'vcf.gz.tbi', 'tsv']])

# generate a cohort level deepvariant vcf or use singleton vcf
include: 'rules/cohort_glnexus.smk'
if 'deepvariant_vcf' in config['cohort_targets']:
    targets.extend([slivar_input, slivar_input + '.tbi'])

# annotate and filter deepvariant vcf
include: 'rules/cohort_slivar.smk'
if 'slivar' in config['cohort_targets']:
    targets.extend([f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.{infix}.{suffix}"
                    for infix in ['slivar', 'slivar.compound-hets']
                    for suffix in ['vcf.gz', 'vcf.gz.tbi', 'tsv']])


ruleorder: split_glnexus_vcf > whatshap_phase > whatshap_bcftools_concat > bcftools_bcf2vcf > bgzip_vcf
localrules: all, md5sum


rule all:
    input: targets + [f"{x}.md5" for x in targets]


rule md5sum:
    input: "{prefix}"
    output: "{prefix}.md5"
    message: "Creating md5 checksum for {input}."
    shell: "md5sum {input} > {output}"
