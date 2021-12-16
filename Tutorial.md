# Tutorial for PacBio Human WGS Workflow

## Table of Contents

- [Workflow Overview](#workflow-overview)
- [Getting Started](#getting-started)
  - [1. Dependencies](#1-dependencies)
  - [2. Prepare Workspace](#2-prepare-workspace)
  - [3. Configuration Files](#3-configuration-files)
  - [4. Run Analysis](#4-run-analysis)
  - [5. Outputs](#5-outputs)
- [Citations](#citations)

---

# Workflow Overview
The purpose of this pipeline is to comprehensively detect and prioritize variants in human genomes using PacBio HiFi reads. It consists of **three [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows** designed to run sequentially, with PacBio HiFi BAMs or FASTQs as the primary input:

1. [process_smrtcells](#1-process_smrtcells)
2. [process_sample](#2-process_sample)
3. [process_cohort](#3-process_cohort)

## 1. process_smrtcells
A single sample will often be sequenced on multiple SMRT Cells, so this workflow aligns HiFi reads from each SMRT Cell to the GRCh38 human reference genome separately. This allows for quality control steps to confirm, for example, that all HiFi reads are from the same sample.

| Tool  | Task  |
| :---   | :---   |
| [pbmm2](https://github.com/PacificBiosciences/pbmm2)  | align HiFi reads (BAMs or FASTQs) to reference (GRCh38 by default) |
| [mosdepth](https://github.com/brentp/mosdepth)        | calculate aligned coverage depth |
| [scripts/extract_read_length_and_qual.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/extract_read_length_and_qual.py) | generate read length and quality statistics |
| [scripts/infer_sex_from_coverage.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/infer_sex_from_coverage.py) | calculate depth ratios (chrX:chrY, chrX:chr2) from mosdepth summary for sample swap detection |
| [scripts/calculate_M2_ratio.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/calculate_M2_ratio.py) | calculate depth ratio (chrM:chr2) from mosdepth summary for sample swap detection |
| [jellyfish](https://github.com/gmarcais/Jellyfish), [scripts/modimer.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/modimer.py) | count kmers in HiFi reads to dump and export modimers for sample swap detection |

## 2. process_sample
The primary goals of this workflow are variant discovery, variant calling, and assembly for each sample. 

| Tool  | Task  |
| :---   | :---   |
| [pbsv](https://github.com/PacificBiosciences/pbsv)              | call structural variants | 
| [DeepVariant](https://github.com/google/deepvariant)            | call small variants |
| [WhatsHap](https://github.com/whatshap/whatshap/)               | phase small variants and generate merged, haplotagged BAM|
| [mosdepth](https://github.com/brentp/mosdepth)        | calculate aligned coverage depth of merged, haplotagged BAM |
| [tandem-genotypes](https://github.com/mcfrith/tandem-genotypes), [scripts/check_tandem_repeat_coverage.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/check_tandem_repeat_coverage.py) | genotype known tandem repeat expansions associated with disease | 
| [hifiasm](https://github.com/chhylp123/hifiasm)                 | assemble reads | 
| [calN50](https://github.com/lh3/calN50)                      | calculate assembly stats |
| [seqtk](https://github.com/lh3/seqtk) | split assembly contigs into 200kb chunks to facilitate visualization of aligned assembly with IGV |
| [minimap2](https://github.com/lh3/minimap2)                     | align assembly to reference (assembly contigs are split into 200kb chunks to facilitate visualization with IGV) |
| [jellyfish](https://github.com/gmarcais/Jellyfish) | merge jellyfish kmer counts for by sample |
| [scripts/check_kmer_consistency.py](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/scripts/check_kmer_consistency.py) | calculate consistency of kmers between sequencing runs for sample swap detection |


## 3. process_cohort
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

---------------

# Getting Started
These workflows are designed to be implemented on a **linux cluster** and require that you follow the installation and configuration procedures described below.

# 1. Dependencies
- [singularity>=3.5.3](#singularity) installed by root
- [conda](#conda)
- [other](#other)
  - lockfile==0.12.2
  - python3
  - snakemake>=5.19
  - mamba (optional, but recommended)

## Singularity
[Singularity](https://sylabs.io/guides/3.6/admin-guide/installation.html#installation-on-linux) **must be installed with root privileges**. If you do not have root privileges and singularity is not installed on your linux cluster, we recommend contacting your HPC administrator for assistance.

## Conda
If conda is not already available to you, we recommend installing Miniconda, a free minimal installer for conda. Download the latest Linux installer for Miniconda3 (64 bit) from [the Miniconda website](https://docs.conda.io/en/latest/miniconda.html#linux-installers). The correct file is named `Miniconda3-latest-Linux-x86_64.sh` where "latest" is replaced by the most recent version number. Then, execute the bash installer with the following command.

`bash Miniconda3-latest-Linux-x86_64.sh`

## Other
Once conda is installed, the easiest way to manage the final dependencies is by creating a conda environment. The following command creates a conda environment named `pacbio-human-wgs` with the final requirements. 

`conda create -n pacbio-human-wgs -c bioconda -c conda-forge lockfile==0.12.2 python=3 snakemake>=5.19 mamba`

This environment must be activated before using the workflow.

`conda activate pacbio-human-wgs`

[Back to top](#TOP)

---------------

# 2. Prepare Workspace
These snakemake workflows require a very specific directory structure in order to function properly. Empty directories that will store input and output files from the analysis were not built into the repo; this allows users to `git pull` the most recent version of the repo without affecting their own data and analysis files. However, this requires that users build the directory structure themselves before using the workflows for the first time.

First, create a clean project directory and move into that directory.
```
mkdir <directory_name>
cd <directory_name>
```

Then, clone the repo and submodules (`--recursive` flag) into a folder named `workflow/`.
```
git clone --recursive https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake.git workflow
```

Create additional empty directories required by the workflow.
```
mkdir -p cluster_logs smrtcells/ready smrtcells/done samples cohorts
```

There are two additional folders (`reference/` and `resources/`) which contain content necessary for these workflows to run. We are working to make these folders available in a public repository, but until then please contact one of the repository contributors or, if applicable, your PacBio representative to request these materials.


After completing these steps, you can visualize the complete directory structure using the `tree -d` command. It is imperative that you preserve this directory structure (and the contents of those directories) to ensure the proper functioning of the workflows. Avoid moving files/directories or changing their names if you intend to run the workflows again in the future.

```
<directory_name>
    ├── cluster_logs
    ├── cohorts
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
    ├── samples
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

[Back to top](#TOP)

---------------

# 3. Configuration Files
Configuration files are written with [`yaml` syntax](https://yaml.org/spec/1.2.2/). 

## Analysis configuration 
The following configuration files require your attention before running the workflows.
- [`cohort.yaml`](#cohortyaml) is a configuration file that **must be created** for the `process_cohort` workflow. It specifies sample names, relationships, and phenotypes that will be included in the analysis. Please create a file called `cohort.yaml` to reflect the specifics of your own analysis. 
- [`config.yaml`](#configyaml) specifies which steps are run in each workflow and contains file paths and version numbers for the docker images used by singularity.
  
### `cohort.yaml`
The `process_cohort` workflow can be configured to run for singletons, trios, or other related/unrelated cohorts of one or more individuals. Multiple singletons and/or cohorts can be listed in the same `cohort.yaml`. Multiple HPO phenotypes ([Human Phenotype Ontology](https://hpo.jax.org/app/) codes) can be listed. The `HP:0000001` phenotype can be used to specify an unknown phenotype for an affected sample. No phenotype needs to be specified if all samples are `unaffecteds`. 

Notes: 1) Choose `cohort_id` and `sample_id` names that make sense to computers! Avoid names with spaces and non-standard symbols. 2) Preserve the white space and structure of `cohort.yaml` entries (e.g. indents, dashes, spaces) to avoid unintended cohort results.
```
# Singleton
- id: <cohort_id>
  phenotypes:
  - HP:0000001
  affecteds:
  - id: singleton-sampleid
    sex: MALE

# Trio
- id: <cohort_id>
  phenotypes:
  - HP:0000001 
  affecteds:
  - id: trio-probandid
    parents:
    - trio-fatherid
    - trio-motherid
    sex: MALE
  unaffecteds:
  - id: trio-fatherid
    sex: MALE
  - id: trio-motherid
    sex: FEMALE
```
#### Example `cohort.yaml`
```
# Trio with affected proband
- id: EpilepsyTrio17  # cohort id
  phenotypes:
    - HP:0001250  # Seizure
    - HP:0001263  # Global developmental delay
  affecteds:
    - id: family17-01  # sample id of the proband
      sex: MALE
      parents: 
      - family17-02  # sample id of the mother
      - family17-03  # sample id of the father
  unaffecteds:
    - id: family17-02  # sample id of the mother
      sex: FEMALE
    - id: family17-03  # sample id of the father
      sex: MALE

# Unrelated cohort with all unaffected
- id: YorubanTrio  # cohort id
  unaffecteds:
    - id: NA19240  # sample id
      sex: FEMALE
    - id: NA19238  # sample id
      sex: FEMALE
    - id: NA19239  # sample id
      sex: MALE

# Singleton affected but unknown phenotype
- id: hg002  # cohort id
  phenotypes:
  - HP:0000001  # example HPO phenotype for "All" phenotype
  affecteds:
  - id: hg002  # sample id
    sex: MALE
```

### `config.yaml`

By default, all steps in a workflow will run when the workflow is launched. If you'd prefer to only run a few steps in one of the workflows, then you can "comment out" the workflow targets in `config.yaml` by adding a `#` symbol at the beginning of the line. Be aware, however, that some steps require the output of other steps in the workflow. For example, the following configuration where `whatshap` has been commented out would cause errors in the `process_sample` and `process_cohort` workflows because `whatshap` output is required by steps like `tandem-genotypes` and `slivar`.

```
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


## Cluster configuration 
The following configuration files **may need to be edited** based on the specifics of your HPC cluster. Additional flags may be necessary for job submission and you may need to talk to your HPC administrator if the default cluster configuration files aren’t working. Warning: unintential or misinformed changes to these files may prevent the workflows from running properly. 
- `*.cluster.yaml` files for each workflow contain example cluster configurations for a **SLURM cluster** with a `compute` queue for general compute and a `ml` queue for GPU
- `*.cluster.sge.yaml` files for each workflow contain example cluster configurations for a **SGE cluster**. Warning: SGE configuration files are not as well maintained as SLURM configuration files in this repo because they are used less frequently.

## Additional configuration files
The following configuration files should not be modified.
- `reference.yaml` contains file paths and names related to the reference genome and resource files used by various steps in the workflows
- Additional conda environment `.yaml` files for individual steps in the workflows are located in `rules/envs/`

[Back to top](#TOP)

---------------

# 4. Run Analysis
The following instructions are specific to a slurm cluster (i.e. `sbatch`). Users of SGE or related job management systems will need to use appropriate job submission execution and flags. 

1. Move into the analysis directory that you created in [Prepare Workspace](#2-prepare-workspace).
```
cd <directory_name>
```

2. Create `cohort.yaml` to reflect the sample names, phenotypes, and relationships of your samples (see [Configuration Files](#3-configuration-files)). An easy way to do this is to copy the `example_cohort.yaml` from `workflow/` your project directory and edit it with your favorite text editor. Remember to choose `cohort_id` and `sample_id` names that make sense to computers! Avoid names with spaces and non-standard symbols, and preserve the existing list structure in each entry (indents, dashes, spaces).

```
cp workflow/example_cohort.yaml cohort.yaml
nano cohort.yaml
```

3. Create a directory for each sample in `smrtcells/ready`. The names of these directories must match the sample IDs specified in `cohort.yaml`.
```
mkdir smrtcells/ready/<sample_id>
```

4. Put PacBio HiFi reads into their respective directories. The easiest way to do this is with a symlink. **Note: unaligned BAM and FASTQ filenames must be identifiable as HiFi reads, i.e. have the following format.** 
   - regex for BAM: `/m\d{5}[Ue]?_\d{6}_\d{6}.(ccs|hifi_reads).bam`
     - example: `m54119U_210108_012126.ccs.bam`
     - example: `m64013e_210917_004210.hifi_reads.bam`
   - regex for FASTQ: `/m\d{5}[Ue]?_\d{6}_\d{6}.fastq.gz`
     - example: `m54119U_210108_012126.fastq.gz`
     - example: `m64013e_210917_004210.fastq.gz`
  
```
ln -s /path/to/HiFi/BAM/or/FASTQ/<hifi_reads_filename> smrtcells/ready/<sample_id>/
```

5. Activate conda environment

```
conda activate pacbio-human-wgs
```

6. Run `process_smrtcells` workflow. This will process all samples located in `smrtcells/ready`. If you have samples in this folder that you don't want to process, move them to `smrtcells/done`.
```
sbatch workflow/process_smrtcells.sh
```

7. When `process_smrtcells.sh` has completed, run the `process_sample.sh` workflow. If you've logged out of your HPC session, make sure to re-activate the conda environment before submitting this job.
```
sbatch workflow/process_sample.sh <sample_id>
```

8. When all of the samples have completed `process_smrtcells` and `process_sample` workflows, you can run the `process_cohort` workflow. The following command will only process the `<cohort_id>` or specified, which must match an entry in `cohort.yaml`. If you've logged out of your HPC session, make sure to re-activate the conda environment before submitting this job.
```
sbatch workflow/process_cohort.sh <cohort_id>
```


---
**Note**: The first time you run any snakemake workflow (`process_smrtcells`, `process_sample`, `process_cohort`), run one sample first so that the conda environments for that workflow are installed only once. After the snakemake workflow starts launching its own sub-jobs, then you can run the rest of your samples.


[Back to top](#TOP)

---------------

# 5. Outputs
The results of `process_smrtcells` and `process_sample` are written to the `samples/` directory. The results of `process_cohort` are written to the `cohorts/` directory. 

## Key outputs in `samples/`
After `process_smrtcells` and `process_sample` have finished, the top level directory structure of each sample in `samples/` should look like the following:

```
$ tree -dL 1 samples/<sample_id>

samples/<sample_id>
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
├── whatshap
└── whatshap_intermediate

13 directories
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


## Key outputs in `cohorts/`
After `process_cohort` has finished, the top level directory structure of each cohort in `cohorts/` should look like the following:
```
$ tree -dL 1 cohorts/<cohort_id>

cohorts/<cohort_id>
├── benchmarks
├── hifiasm   # only produced if trios present in cohort
├── logs
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

---------------

# Citations
The following resources were used for variant prioritization and should be cited if this workflow generates results for publication. 

| Resource | Citation |
| :--- | :--- |
| Evan Eichler Lab |Audano, Peter A., et al. "Characterizing the major structural variant alleles of the human genome." Cell 176.3 (2019): 663-675. |
| deCODE genetics | Beyter, Doruk, et al. "Long-read sequencing of 3,622 Icelanders provides insight into the role of structural variants in human diseases and other traits." Nature Genetics 53.6 (2021): 779-786. |
| gnomAD-SV | Collins, Ryan L., et al. "An open resource of structural variation for medical and population genetics." BioRxiv (2019): 578674. |
| Human Pangenome Reference Consortium | https://github.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0 |
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
| pbmm2 | https://github.com/PacificBiosciences/pbmm2 |
| pbsv | https://github.com/pacificbiosciences/pbsv |
| slivar | Pedersen, Brent S., et al. "Effective variant filtering and expected candidate variant yield in studies of rare human disease." NPJ Genomic Medicine 6.1 (2021): 1-8. |
| Snakemake | Mölder, Felix, et al. "Sustainable data analysis with Snakemake." F1000Research 10 (2021). |
| svpack | https://github.com/PacificBiosciences/svpack |
| tandem-genotypes | Mitsuhashi, Satomi, et al. "Tandem-genotypes: robust detection of tandem repeat expansions from long DNA reads." Genome biology 20.1 (2019): 1-17. |
| WhatsHap | Martin, Marcel, et al. "WhatsHap: fast and accurate read-based phasing." BioRxiv (2016): 085050. |


[Back to top](#TOP)
