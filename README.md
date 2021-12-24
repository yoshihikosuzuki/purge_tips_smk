# purge_tips_smk

Snakemake version of [purge_tips](https://github.com/yoshihikosuzuki/purge_tips).

## Dependency

- Snakemake
- samtools
- Hifiasm
- gfatools
- Meryl
- Winnowmap

We strongly recommend installing these software 1) as Environment Module/Lmod modules or 2) to some location accessible via `PATH`; if so, then you can easy specify the dependencies in the config file of purge_tips.

## How to use/install

Essentially you can run purge_tips simply by executing the `purge_tips` file (which is actually a bash script invoking the `snakemake` command) in the `workflow/` directory in this repository:

```bash
$ git clone https://github.com/yoshihikosuzuki/purge_tips_smk
$ cd purge_tips_smk
$ purge_tips [snakemake options]
```

To make the command available at any directory, add the root path to the `PATH` environment variable:

```bash
$ export PATH=/path/to/purge_tips_smk:$PATH
```

## How to specify inputs/parameters/other configs

All configurations are spcified via a config file, whose template is the `config.yaml` file in this repository:

```bash
################################################################################
#                        Environment-specific settings
################################################################################

# Commands to load modulefiles, etc.
shell_prefix: >-
  shopt -s expand_aliases;
  . /apps/free/lmod/lmod/init/bash;
  module use /apps/.bioinfo-ugrp-modulefiles81;

# List of dependent environment modules
modules:
  seqkit: "Other/seqkit/2.0.0"
  samtools: "samtools/1.12"
  hifiasm: "Other/hifiasm/0.16.1"
  gfatools: "Other/gfatools/0.5"
  winnowmap: "Other/winnowmap/2.03"

################################################################################
#                         Dataset-specific settings
################################################################################

hifi_fastq: "your-hifi-fastq-file"
error_thres: 10
hifiasm_options: ""
```
