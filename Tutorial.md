# Tutorial for PacBio Human WGS Workflow

## Table of Contents

- [Workflow Overview](#workflow-overview)
- [Getting Started](#getting-started)
  - [1. Dependencies](#1-dependencies)
  - [2. Prepare Workspace](#2-prepare-workspace)
  - [3. Analysis Configuration](#3-analysis-configuration)
  - [4. Cluster Configuration](#4-cluster-configuration)
  - [5. Run Analysis](#5-run-analysis)
- [Outputs](#outputs)
- [Citations](#citations)

---

## Workflow Overview

The purpose of this pipeline is to comprehensively detect and prioritize variants in human genomes using PacBio HiFi reads. It consists of **three [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows** designed to run sequentially, with PacBio HiFi BAMs or FASTQs as the primary input:

1. [process_smrtcells](#1-process_smrtcells)
2. [process_sample](#2-process_sample)
3. [process_cohort](#3-process_cohort)

### process_smrtcells

A single sample will often be sequenced on multiple SMRT Cells, so this workflow aligns HiFi reads from each SMRT Cell to the GRCh38 human reference genome separately. This allows for quality control steps to confirm, for example, that all HiFi reads are from the same sample.

| Tool  | Task  |
| :---   | :---   |
| [pbmm2](https://github.com/PacificBiosciences/pbmm2)  | align HiFi reads (BAMs or FASTQs) to reference (GRCh38 by default) |
| [mosdepth](https://github.com/brentp/mosdepth)        | calculate aligned coverage depth |
| [scripts/extract_read_length_and_qual.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/extract_read_length_and_qual.py) | generate read length and quality statistics |
| [scripts/infer_sex_from_coverage.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/infer_sex_from_coverage.py) | calculate depth ratios (chrX:chrY, chrX:chr2) from mosdepth summary for sample swap detection |
| [scripts/calculate_M2_ratio.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/calculate_M2_ratio.py) | calculate depth ratio (chrM:chr2) from mosdepth summary for sample swap detection |
| [jellyfish](https://github.com/gmarcais/Jellyfish), [scripts/modimer.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/modimer.py) | count kmers in HiFi reads to dump and export modimers for sample swap detection |

### process_sample

The primary goals of this workflow are variant discovery, variant calling, and assembly for each sample.

| Tool  | Task  |
| :---   | :---   |
| [pbsv](https://github.com/PacificBiosciences/pbsv)              | call structural variants |
| [DeepVariant](https://github.com/google/deepvariant)            | call small variants |
| [pb-cpg-tools](https://github.com/PacificBiosciences/pb-CpG-tools) | obtain list of CpG/5mC sites and modification probabilities |
| [WhatsHap](https://github.com/whatshap/whatshap/)               | phase small variants and generate merged, haplotagged BAM|
| [BCFtools RoH](https://samtools.github.io/bcftools/howtos/roh-calling.html)| detect regions of autozygosity in merged, haplotagged BAM using a hidden Markov model |
| [mosdepth](https://github.com/brentp/mosdepth)        | calculate aligned coverage depth of merged, haplotagged BAM |
| [tandem-genotypes](https://github.com/mcfrith/tandem-genotypes), [scripts/check_tandem_repeat_coverage.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/check_tandem_repeat_coverage.py) | genotype known tandem repeat expansions associated with disease |
| [trgt](https://github.com/PacificBiosciences/trgt) | genotype tandem repeats |
| [hifiasm](https://github.com/chhylp123/hifiasm)                 | assemble reads |
| [calN50](https://github.com/lh3/calN50)                      | calculate assembly stats |
| [seqtk](https://github.com/lh3/seqtk) | split assembly contigs into 200kb chunks to facilitate visualization of aligned assembly with IGV |
| [minimap2](https://github.com/lh3/minimap2)                     | align assembly to reference (assembly contigs are split into 200kb chunks to facilitate visualization with IGV) |
| [jellyfish](https://github.com/gmarcais/Jellyfish) | merge jellyfish kmer counts for by sample |
| [scripts/check_kmer_consistency.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/check_kmer_consistency.py) | calculate consistency of kmers between sequencing runs for sample swap detection |

### process_cohort

In this workflow, variants are prioritized, annotated, and filtered to aid in the process of finding candidate rare variants with functional consequence. It can be run with a single sample (singleton) or as a multi-sample cohort (e.g., trio, quad) depending on the configuration specified in `cohort.yaml` ([details below](#configuration-files)).

| Tool  | Task  |
| :---   | :---   |
| [pbsv](https://github.com/PacificBiosciences/pbsv) | joint call structural variants (**for multi-sample cohorts**) |
| [GLnexus](https://github.com/dnanexus-rnd/GLnexus) | joint call small variants (**for multi-sample cohorts**) |
| [slivar](https://github.com/brentp/slivar) | annotate small variant calls with population frequency from [gnomAD](https://gnomad.broadinstitute.org) and [HPRC](https://humanpangenome.org) variant databases |
| [bcftools](https://github.com/samtools/bcftools) | annotate small variant calls with functional consequence based on position and reference genome annotation |
| [slivar](https://github.com/brentp/slivar) | filter small variant calls according to population frequency and inheritance patterns |
| [slivar](https://github.com/brentp/slivar), [scripts/add_comphet_phase.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/add_comphet_phase.py) | detect possible compound heterozygotes, and filter to remove cis-combinations |
| [scripts/calculate_phrank.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/calculate_phrank.py) | assign a phenotype rank ([Phrank](https://www.nature.com/articles/s41436-018-0072-y)) score to variants |
| [svpack](https://github.com/PacificBiosciences/svpack) | annotate structural variant calls with population frequency from [gnomAD-SV](https://gnomad.broadinstitute.org/), [HPRC](https://humanpangenome.org), and additional published datasets (see [resources](#citations)) |
| [svpack](https://github.com/PacificBiosciences/svpack) | filter structural variant calls based on population frequency  |
| [svpack](https://github.com/PacificBiosciences/svpack) | annotate structural variant calls with functional consequence based on position and reference genome annotation |
| [slivar](https://github.com/brentp/slivar) | convert filtered and annotated variants from slivar and svpack into TSVs for easy browsing |
| [yak](https://github.com/lh3/yak) | create k-mer database for each parent (only if trios included in cohort) |
| [hifiasm](https://github.com/chhylp123/hifiasm)                 | assemble reads with trio binning (only if trios included in cohort) |
| [calN50](https://github.com/lh3/calN50)                      | calculate assembly stats  (only if trios included in cohort) |
| [seqtk](https://github.com/lh3/seqtk) | split assembly contigs into 200kb chunks to facilitate visualization of aligned assembly with IGV (only if trios included in cohort) |
| [minimap2](https://github.com/lh3/minimap2)                     | align assembly to reference (only if trios included in cohort) |

[Back to top](#TOP)

---

## Getting Started

These workflows are designed to be implemented on a **linux cluster** and require that you follow the installation and configuration procedures described below.

### 1. Dependencies

- [singularity>=3.5.3](#singularity) installed by root
  - [Singularity](https://sylabs.io/guides/3.6/admin-guide/installation.html#installation-on-linux) **must be installed with root privileges**. If you do not have root privileges and singularity is not installed on your linux cluster, we recommend contacting your HPC administrator for assistance.
- [conda](#conda)
  - If conda is not already available to you, we recommend installing Miniconda, a free minimal installer for conda. Download the latest Linux installer for Miniconda3 (64 bit) from [the Miniconda website](https://docs.conda.io/en/latest/miniconda.html#linux-installers).

[Back to top](#TOP)

---

### 2. Prepare Workspace

These snakemake workflows require a very specific directory structure in order to function properly. Empty directories that will store input and output files from the analysis were not built into the repo; this allows users to `git pull` the most recent version of the repo without affecting their own data and analysis files. However, this requires that users build the directory structure themselves before using the workflows for the first time.

```text
# create a clean project directory and move into that directory
mkdir <directory_name>
cd <directory_name>

# clone the repo and submodules (--recursive flag) into a folder named workflow/
git clone --recursive https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake.git workflow

# create additional empty directories required by the workflow
mkdir -p cluster_logs smrtcells/ready smrtcells/done samples cohorts
```

There are two additional folders (`reference/` and `resources/`) which contain content necessary for these workflows to run. We are working to make these folders available in a public repository, but until then please contact one of the repository contributors or, if applicable, your PacBio representative to request these materials.

After completing these steps, you can visualize the complete directory structure using the `tree -d` command. It is imperative that you preserve this directory structure (and the contents of those directories) to ensure the proper functioning of the workflows. Avoid moving files/directories or changing their names if you intend to run the workflows again in the future.

```text
<directory_name>
    ├── cluster_logs
    ├── cohorts  # created during process_cohort
    ├── reference
    │   └── annotation
    ├── resources
    │   ├── decode
    │   ├── eee
    │   ├── gnomad
    │   ├── gnomadsv
    │   ├── hpo
    │   ├── hprc
    │   ├── jellyfish
    │   ├── slivar
    │   └── tandem-genotypes
    ├── samples  # created during process_smrtcells
    ├── smrtcells
    │   ├── done
    │   └── ready
    └── workflow
        ├── rules
        │   └── envs
        └── scripts
            ├── calN50
            └── svpack
```

Now you can add your input files (unaligned PacBio HiFi reads) to the correct location.

```text
# create a directory for each sample in smrtcells/ready.
mkdir smrtcells/ready/<sample_id>


# add one or more unaligned PacBio HiFi reads files in the correct sample directory with a symlink 
ln -s /path/to/HiFi/BAM/or/FASTQ/<hifi_reads_filename> smrtcells/ready/<sample_id>/
```

> **WARNING: unaligned BAM and FASTQ filenames must be identifiable as HiFi reads, i.e. have the following format.**
>
> - regex for BAM: `/m\d{5}[Ue]?_\d{6}_\d{6}.(ccs|hifi_reads).bam`
>   - example: `m54119U_210108_012126.ccs.bam`
>   - example: `m64013e_210917_004210.hifi_reads.bam`
> - regex for FASTQ: `/m\d{5}[Ue]?_\d{6}_\d{6}.fastq.gz`
>   - example: `m54119U_210108_012126.fastq.gz`
>   - example: `m64013e_210917_004210.fastq.gz`

[Back to top](#TOP)

---

### 3. Analysis Configuration

Configuration files are written with [`yaml` syntax](https://yaml.org/spec/1.2.2/).

#### Create `cohort.yaml`

The `process_cohort` workflow requires you to specify sample names, relationships, and phenotypes in `cohort.yaml`. Cohorts can consist of singletons, trios, or other related/unrelated groups of one or more individuals. Multiple singletons and/or cohorts can be listed in the same `cohort.yaml` and run independently. Multiple HPO phenotypes ([Human Phenotype Ontology](https://hpo.jax.org/app/) codes) can be listed for each sample. The `HP:0000001` term can be used to specify an unknown phenotype for an affected sample. No phenotype needs to be specified if all samples are `unaffecteds`.

The layout of `cohort.yaml` is shown below and an example `cohort.yaml` file can be found [here](example_cohort.yaml).

```text
# Singleton
- id: <cohort_id>
  phenotypes:
  - HP:0000001
  affecteds:
  - id: <sample_id>
    sex: MALE

# Trio
- id: <cohort_id>
  phenotypes:
  - HP:0000001 
  affecteds:
  - id: <child_id>
    parents:
    - <father_id>
    - <mother_id>
    sex: MALE
  unaffecteds:
  - id: <father_id>
    sex: MALE
  - id: <mother_id>
    sex: FEMALE
```

> **WARNING**:
>
> - Choose `cohort_id` and `sample_id` names that make sense to computers! Avoid names with spaces and non-standard symbols.
> - Each cohort must have a unique `cohort_id`.
> - Preserve the white space and structure of `cohort.yaml` entries (e.g. indents, dashes, spaces) to avoid unintended cohort results.

#### Modify `config.yaml` if necessary

The `config.yaml` file specifies which steps are run in each workflow and contains various parameters, file paths, and version numbers for the docker images used by singularity. By default, all steps in a workflow will run when the workflow is launched. If you'd prefer to only run a few steps in one of the workflows, then you can "comment out" the workflow targets in `config.yaml` by adding a `#` symbol at the beginning of the line. Be aware, however, that some steps require the output of other steps in the workflow. For example, the following configuration where `whatshap` has been commented out would cause errors in the `process_sample` and `process_cohort` workflows because `whatshap` output is required by steps like `tandem-genotypes` and `slivar`.

```text
smrtcells_targets:
  - alignment
  - stats  # req: alignment
  - coverage  # req: alignment
  - coverage_qc  # req: alignment
  - kmers

sample_targets:
  - pbsv_vcf  # req: alignment in config['smrtcells_targets']
  - deepvariant  # req: alignment in config['smrtcells_targets']
#  - whatshap  # req: deepvariant
  - coverage  # req: whatshap
  - kmers  # req: kmers in config['smrtcells_targets']
  - assembly
  - tandem-genotypes  # req: whatshap

cohort_targets:
  - pbsv_vcf  # req: pbsv_vcf in config['sample_targets']
  - svpack  # req: pbsv_vcf in config['sample_targets']
  - deepvariant_vcf  # req: deepvariant, whatshap in config['sample_targets']
  - slivar  # req: deepvariant, whatshap in config['sample_targets']
  - trio_assembly
```

#### Additional configuration files

The following configuration files should not be modified unless you're confident about what you're doing.

- `reference.yaml` contains file paths and names related to the reference genome and resource files used by various steps in the workflows
- Additional conda environment `.yaml` files for individual steps in the workflows are located in `rules/envs/`

[Back to top](#TOP)

---

### 4. Cluster configuration

We have provided sample submission scripts for three different schedulers (SLURM, SGE, LSF) as well as a local execution option. The following files **may need to be edited** based on the specifics of your HPC cluster. Additional flags may be necessary for job submission and you may need to talk to your HPC administrator if the example scripts aren’t working.

- `process_smrtcells.(slurm|sge|lsf|local).sh`
- `process_sample.(slurm|sge|lsf|local).sh`
- `process_cohort.(slurm|sge|lsf|local).sh`
- `profiles/(slurm|sge|lsf|local)/config.yaml` # snakemake configuration [profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
- `rules/sample_deepvariant.smk` # GPU and singularity constraints specified within rules
- `variables.env` # environment variables such as TMPDIR, cluster partition, and cluster account

> **WARNING**:
>
> - We recommend at least 80 cores and 1TB RAM for local execution.  Local execution will use all available cores.
> - Job scripts and configuration for SGE and LSF schedulers are included as a courtesy, but are not regularly tested or maintained.
> - Unintential or misinformed changes to these files may prevent the workflows from running properly.
> - If a cluster "account" is required for proper billing, change the account value from `100humans` in `variables.env` and in the relevant `*.(slurm|sge|lsf|local).sh` files.

[Back to top](#TOP)

---

### 5. Run Analysis

First, create a conda environment with the required packages. This environment only needs to be created once.

```text
# create conda environment
conda install mamba -n base -c conda-forge
conda activate base
mamba create -c conda-forge -c bioconda -n pb-human-wgs snakemake=6.15.3 tabulate=0.8.10 pysam=0.16.0.1 python=3
```

To run the workflows on a cluster that uses the Slurm job scheduler, use the following commands. Users of SGE, LSF, or related job management systems will need to use appropriate job submission execution and flags.

```text
# activate conda environment
# do this every time you want to run any part of workflow
conda activate pb-human-wgs

# confirm that you're in the analysis directory you created (parent directory of workflow)

# run process_smrtcells on all samples in smrtcells/ready
sbatch workflow/process_smrtcells.slurm.sh

# process_smrtcells must finish before launching next step!
# run process_sample on a single sample
sbatch workflow/process_sample.slurm.sh <sample_id>

# process_smrtcells & process_sample must finish for all samples in cohort before next step!
# run process_cohort on a single cohort described in cohort.yaml
sbatch workflow/process_cohort.slurm.sh <cohort_id>
```

For local execution, use the following commands.

```text
# activate conda environment
# do this every time you want to run any part of workflow
conda activate pb-human-wgs

# confirm that you're in the analysis directory you created (parent directory of workflow)

# run process_smrtcells on all samples in smrtcells/ready
bash workflow/process_smrtcells.local.sh

# process_smrtcells must finish before launching next step!
# run process_sample on a single sample
bash workflow/process_sample.local.sh <sample_id>

# process_smrtcells & process_sample must finish for all samples in cohort before next step!
# run process_cohort on a single cohort described in cohort.yaml
bash workflow/process_cohort.local.sh <cohort_id>
```

> **WARNING**:
>
> - The conda environment must activated before launching any of the workflows, so if you've logged out of your HPC session, make sure to re-activate the conda environment before submitting the next job.
> - The first time you run any snakemake workflow (`process_smrtcells`, `process_sample`, `process_cohort`), run one sample first so that the conda environments for that workflow are installed only once. After the snakemake workflow starts launching its own sub-jobs, then you can run the rest of your samples.
> - The `process_smrtcells` workflow will process all samples located in `smrtcells/ready`. If you have samples in this folder that you don't want to process, move them to `smrtcells/done`.
>

[Back to top](#TOP)

---

## Outputs

The results of `process_smrtcells` and `process_sample` are written to the `samples/` directory. The results of `process_cohort` are written to the `cohorts/` directory.

### Key outputs in `samples/`

After `process_smrtcells` and `process_sample` have finished, the top level directory structure of each sample in `samples/` should look like the following:

```text
$ tree -dL 1 samples/<sample_id>

samples/<sample_id>
├── 5mc_cpg_pileup  # only if HiFi reads are provided in BAM format and contain methylation tags (Ml)
├── aligned
├── benchmarks
├── deepvariant
├── deepvariant_intermediate
├── hifiasm
├── jellyfish
├── logs
├── mosdepth
├── pbsv
├── smrtcell_stats
├── tandem-genotypes
├── trgt
├── whatshap
└── whatshap_intermediate

15 directories
```

The following are some of the key output files from these workflows.

- **SMRT Cell Statistics**
  - `smrtcell_stats/*.read_length_and_quality.tsv`
  - One summary file per SMRT Cell with 3 columns:
    - [1]read ID
    - [2]read length (bp)
    - [3]read quality
- **Coverage**
  - `mosdepth/*.GRCh38.mosdepth.summary.txt`
    - Coverage summary statistics
    - One summary file per SMRT Cell and one combined summary file for the sample with 6 columns:
      - [1]chrom: chromosome
      - [2]length: chromosome length
      - [3]bases: number of bases mapped to chromosome
      - [4]mean: average number of reads mapped to each base in chromosome
      - [5]min: minimum number of reads mapped to each base in chromosome
      - [6]max: maximum number of reads mapped to each base in chromosome
  - `mosdepth/*.GRCh38.mosdepth.inferred_sex.txt`
    - Ratios of chromX:chromY coverage and chromX:chrom2 coverage and sex inferred from each
    - Important for determining if a sample swap has occurred in a SMRT Cell or library!
    - One file per SMRT Cell and one combined file for the sample
  - `mosdepth/*.GRCh38.mosdepth.M2_ratio.txt`
    - Contains a single value, the ratio of mtDNA coverage to chrom 2 coverage (chrM:chr2)
    - One file per SMRT Cell
  - `mosdepth/*.GRCh38.gc_coverage.summary.txt`
    - Coverage in regions with X% GC content
    - One file per SMRT Cell and one combined file for the sample
- **Aligned Reads (BAMs)**
  - `whatshap/*.GRCh38.deepvariant.haplotagged.bam` and `.bai`
    - Reads from all SMRT Cells combined, aligned to the reference genome, and tagged with haplotype (HP) determined through small variant phasing + associated index file (.bai)
  - `aligned/*.GRCh38.bam` and `.bai`
    - Reads aligned to the reference genome + associated index file (.bai)
    - One per SMRT Cell
- **Variant Calls**
  - `whatshap/*.GRCh38.deepvariant.phased.vcf.gz` and `.tbi`
    - Phased small variants (SNPs and INDELs < 50bp) + associated index file (.tbi)
  - `deepvariant/*.GRCh38.deepvariant.vcf.stats.txt`
    - Summary statistics for small variants
  - `pbsv/*.GRCh38.pbsv.vcf.gz` and `tbi`
    - Structural variants + associated index file
  - `tandem-genotypes/*.tandem-genotypes.txt`
    - Genotypes for a list of known disease-causing tandem repeat variants
  - `tandem-genotypes/*.tandem-genotypes.absolute.txt`
    - Same as above, but with repeat counts adjusted by adding estimated number of reads in the reference
  - `tandem-genotypes/*.tandem-genotypes.dropouts.txt`
    - Regions with insufficient coverage to genotype are listed here
  - `trgt/*.trgt.vcf.gz` and `tbi`
    - Genotypes for ~170k STR loci
  - `trgt/*.trgt.spanning.bam`
    - Fragments of HiFi reads spanning ~170k STR loci
  - `trgt/*.trgt.dropouts.txt`
    - Regions with insufficient coverage to genotype are listed here
- **Assembly**
  - `hifiasm/*.asm.bp.hap1.p_ctg.fasta.gz`
  - `hifiasm/*.asm.bp.hap2.p_ctg.fasta.gz`
    - Genome assembly split into pseudohaplotypes (hap1 and hap2)
  - `hifiasm/*.asm.bp.hap1.p_ctg.noseq.gfa`
  - `hifiasm/*.asm.bp.hap2.p_ctg.noseq.gfa`
    - Genome assembly graphs for each pseudohaplotype
  - `hifiasm/*.asm.bp.hap1.p_ctg.fasta.stats.txt`
  - `hifiasm/*.asm.bp.hap2.p_ctg.fasta.stats.txt`
    - Summary stats for each pseudohaplotype
  - `hifiasm/*.asm.GRCh38.bam` and `.bai`
    - Genome assembly aligned to the reference + associated index file (.bai)
  - `hifiasm/*.asm.GRCh38.htsbox.vcf.gz` and `.tbi`
    - Variants called after aligned assembly to the reference + associated index file (.tbi)
- **CpG/5mC Scores**
  - `5mc_cpg_pileup/*.denovo.bed` and `5mc_cpg_pileup/*.denovo.bw`
    - Bed and bigwig files for the complete read set and each haplotype showing methylation probabilities for CG sites with a minimum coverage of 4 (default)
  - `5mc_cpg_pileup/*.denovo.mincov10.bed` and `5mc_cpg_pileup/*.denovo.mincov10.bw`
    - Bed and bigwig files for the complete read set and each haplotype showing methylation probabilities for CG sites with a minimum coverage of 10

### Key outputs in `cohorts/`

After `process_cohort` has finished, the top level directory structure of each cohort in `cohorts/` should look like the following:

```text
$ tree -dL 1 cohorts/<cohort_id>

cohorts/<cohort_id>
├── benchmarks
├── glnexus  # only produced if cohort consists of >1 sample
├── hifiasm   # only produced if trios present in cohort
├── logs
├── pbsv  # only produced if cohort consists of >1 sample
├── slivar
└── svpack

5 directories
```

The following are some of the key output files from these workflows. The haplotype-resolved assembly is only produced when a cohort includes one or more trios (child and both parents).

- **Filtered & Annotated Small Variants**
  - `slivar/*.GRCh38.deepvariant.phased.slivar.vcf.gz`
    - Small variant calls that are filtered based on population frequency and annotated with cohort information, population frequency, gene, functional impact, etc.
    - A corresponding `.tsv` is also available in this directory, which can be opened in Excel or other spreadsheet softwares
  - `slivar/*GRCh38.deepvariant.phased.slivar.compound-hets.vcf.gz`
    - Compound heterozygotes annotated with cohort information, population frequency, gene, functional impact, etc.
    - A corresponding `.tsv` is also available in this directory, which can be opened in Excel or other spreadsheet softwares
- **Filtered & Annotated Structural Variants**
  - `svpack/*.GRCh38.pbsv.svpack.vcf.gz` and `.tbi`
    - Structural variant calls that are filtered based on population frequency and annotated with cohort information, population frequency, gene, functional impact, etc.
    - A corresponding `.tsv` is also available in this directory, which can be opened in Excel or other spreadsheet softwares
- **Haplotype-Resolved Assembly**
  - `hifiasm/*.asm.dip.hap1.p_ctg.fasta.gz`
  - `hifiasm/*.asm.dip.hap2.p_ctg.fasta.gz`
    - Genome assembly split into haplotypes (hap1 and hap2)
  - `hifiasm/*.asm.dip.hap1.p_ctg.noseq.gfa`
  - `hifiasm/*.asm.dip.hap2.p_ctg.noseq.gfa`
    - Genome assembly graphs for each haplotype
  - `hifiasm/*.asm.bp.hap1.p_ctg.fasta.stats.txt`
  - `hifiasm/*.asm.bp.hap2.p_ctg.fasta.stats.txt`
    - Summary stats for each haplotype
  - `hifiasm/*.asm.GRCh38.bam` and `.bai`
    - Genome assembly aligned to the reference + associated index file (.bai)

[Back to top](#TOP)

---

## Troubleshooting and FAQ

This section includes problems frequently encountered by users of this pipeline. Please file a repo issue or contact one of the repo contributors if the following troubleshooting tips don't address your concerns.

**Problem:** Workflow won't run and gives error `snakemake: command not found`  
**Solution:** Make sure the conda environment is installed and activated before trying to run the workflow. Instructions [here](#5-run-analysis).

**Problem:** Workflow won't run and gives error `lockfile $LOCKFILE already exists. Remove lockfile and try again.`  
**Solution:** If the workflows are still running or a job failed, the lockfile may not have been properly removed. Either wait for the workflow to finish or, if the job failed, manually remove the lockfile. This error can also be caused if the input file folder doesn't exist. For example, if you try to run `process_sample` without first running `process_smrtcells`.

**Problem:** Trio assembly has one haplotype that is significantly larger than the other  
**Solution:** It's possible that parental HiFi coverage for at least one parent/haplotype is insufficient for the trio binning step in hifiasm. Consider sequencing parents to greater depth.

**Problem:** The `process_smrtcells` workflow starts, but no jobs are executed and it says `uBAMs available for samples: [] FASTQs available for samples: []`  
**Solution:**  Make sure you've provided input files in `smrtcells/ready/<sample_id>` The folder <sample_id> must be created with `mkdir` (not symlinked) although files inside this folder can be symlinked.

**Problem:** I don't have access to GPUs  
**Solution:**  Make sure the line `cpu_only: True` is in `workflow/config.yaml`. Additional changes may be required in `workflow/variables.env` depending on your job scheduler. You might also need to reduce `max-threads` in `workflow/profiles/<profile>/config.yaml` if you don't have access to a 256-core machine.

**Problem:** No space left on device  
**Solution:** You may need to clean out `/tmp` on the host that produced this error. If this issue persists, change the TMPDIR variable in `workflow/variables.env` to a directory that has sufficient space for temporary files.

**Question:** Can I use the CHM13 (T2T) reference genome?  
**Answer:** We don't officially support the chm13v2.0 reference and we haven't developed all of the tertiary analysis resources (variant frequency databases, segdup/repeat/oddregion bed files, etc) to accompany this reference. Even though we don't support it, here is some rough code to get you started:

1. Download the recommended reference fasta [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz) and index with `samtools faidx chm13v2.0_maskedY_rCRS.fa.gz`.
2. Create a chromosome length file from index for phasing with `cut -f1,2 chm13v2.0_maskedY_rCRS.fa.fai > chm13v2.0_maskedY_rCRS.chr_lengths.txt`. Drop the CHM13 reference, index, and chr_lengths file in the `reference/` folder, and update `fasta`, `index`, and `chr_lengths` paths in `workflow/reference.yaml` to match.
3. Download the correct tandem repeats bed file for structural variant calling from [here](https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_chm13v2.0_maskedY_rCRS.trf.bed), drop in the `reference/` folder, and update the `tr_bed` path in `workflow/reference.yaml` to match.
4. In `config.yaml`, disable `kmers` and `tandem-genotypes` under `sample_targets` and `svpack` and `slivar` under `cohort_targets`.
5. See the [CHM13 repo](https://github.com/marbl/CHM13#downloads), for additional resources that might support your analysis, including liftoff annotations, ensembl, refseq, ClinVar, and dbSNP.

[Back to top](#TOP)

---

## Citations

The following resources were used for variant prioritization and should be cited if this workflow generates results for publication.

| Resource | Citation |
| :--- | :--- |
| Evan Eichler Lab |Audano, Peter A., et al. "Characterizing the major structural variant alleles of the human genome." Cell 176.3 (2019): 663-675. |
| deCODE genetics | Beyter, Doruk, et al. "Long-read sequencing of 3,622 Icelanders provides insight into the role of structural variants in human diseases and other traits." Nature Genetics 53.6 (2021): 779-786. |
| gnomAD-SV | Collins, Ryan L., et al. "An open resource of structural variation for medical and population genetics." BioRxiv (2019): 578674. |
| Human Pangenome Reference Consortium | [HPRC Year1 Data Freeze](https://github.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0) |
| Phrank Scores | Jagadeesh, Karthik A., et al. "Phrank measures phenotype sets similarity to greatly improve Mendelian diagnostic disease prioritization." Genetics in Medicine 21.2 (2019): 464-470. |
| gnomAD | Karczewski, K.J., Francioli, L.C., Tiao, G. et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 581, 434–443 (2020). doi:10.1038/s41586-020-2308-7 |
| Human Phenotype Ontology | Köhler, Sebastian, et al. "The human phenotype ontology in 2021." Nucleic acids research 49.D1 (2021): D1207-D1217. |
| pathogenic tandem repeat expansions | Tang, Haibao et al. “Profiling of Short-Tandem-Repeat Disease Alleles in 12,632 Human Whole Genomes.” American journal of human genetics vol. 101,5 (2017): 700-715. doi:10.1016/j.ajhg.2017.09.013 |

The following tools and methods should also be cited if this workflow generates results for publication.
| Tool | Citation |
| :--- | :--- |
| BEDTools | Quinlan, Aaron R., and Ira M. Hall. "BEDTools: a flexible suite of utilities for comparing genomic features." Bioinformatics 26.6 (2010): 841-842. |
| BCFtools, SAMtools | Danecek P, Bonfield JK, et al. Twelve years of SAMtools and BCFtools. Gigascience (2021) 10(2):giab008 |
| DeepVariant | Poplin, Ryan, et al. "A universal SNP and small-indel variant caller using deep neural networks." Nature biotechnology 36.10 (2018): 983-987. |
| GLnexus | Yun, Taedong, et al. "Accurate, scalable cohort variant calls using DeepVariant and GLnexus." Bioinformatics 36.24 (2020): 5582-5589. |
| hifiasm | Cheng, Haoyu, et al. "Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm." Nature Methods 18.2 (2021): 170-175. |
| HTSlib | Bonfield, James K., et al. "HTSlib: C library for reading/writing high-throughput sequencing data." Gigascience 10.2 (2021): giab007. |
| Jellyfish | Marçais, Guillaume, and Carl Kingsford. "A fast, lock-free approach for efficient parallel counting of occurrences of k-mers." Bioinformatics 27.6 (2011): 764-770. |
| LAST | Kiełbasa, Szymon M., et al. "Adaptive seeds tame genomic sequence comparison." Genome research 21.3 (2011): 487-493. |
| minimap2 | Li, Heng. "Minimap2: pairwise alignment for nucleotide sequences." Bioinformatics 34.18 (2018): 3094-3100. |
| mosdepth | Pedersen, Brent S., and Aaron R. Quinlan. "Mosdepth: quick coverage calculation for genomes and exomes." Bioinformatics 34.5 (2018): 867-868. |
| pbmm2 | [pbmm2 GitHub repo](https://github.com/PacificBiosciences/pbmm2) |
| pbsv | [pbsv GitHub repo](https://github.com/pacificbiosciences/pbsv) |
| slivar | Pedersen, Brent S., et al. "Effective variant filtering and expected candidate variant yield in studies of rare human disease." NPJ Genomic Medicine 6.1 (2021): 1-8. |
| Snakemake | Mölder, Felix, et al. "Sustainable data analysis with Snakemake." F1000Research 10 (2021). |
| svpack | [svpack GitHub repo](https://github.com/PacificBiosciences/svpack) |
| tandem-genotypes | Mitsuhashi, Satomi, et al. "Tandem-genotypes: robust detection of tandem repeat expansions from long DNA reads." Genome biology 20.1 (2019): 1-17. |
| trgt | [trgt GitHub repo](https://github.com/PacificBiosciences/trgt) |
| WhatsHap | Martin, Marcel, et al. "WhatsHap: fast and accurate read-based phasing." BioRxiv (2016): 085050. |

[Back to top](#TOP)
