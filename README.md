# purge_tips

Given an assembly (as a GFA file) and reads from which the assembly was generated, purge_tips removes reads mapped to low-coverage contigs that are deemed contamination, etc.
Such low-coverage contigs often exist as "tips" (i.e. short branches) in an assembly graph, and we named this tool after [purge_dups](https://github.com/dfguan/purge_dups), one of whose purposes is somewhat similar to ours.
By simply re-running an assembler with the filtered reads, it is expected that one can possibly improve the N50 and the quality of the assembly.

One use case would be when the total size of the assembled primary unitigs (which should represent all the sequences from which the reads were generated, without merging different haplotypes) is too large compared to the known genome size times the ploidy of the genome, i.e. the maximum size of an ideal genome assembly of a single individual.

Currently we assume only HiFi reads as input and assume that there exist some environment modules with specific names (hard-coded in `purge_tips`).

**Disclaimer: This is experimental at this point, and we have not confirmed if it is accurate in any sense.**

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
$ purge_tips [snakemake options]
```

To make the command available at any directory, add the path to `workflow/` to the `PATH` environment variable:

```bash
$ export PATH=/path/to/workflow/:$PATH
```

## How to specify inputs/parameters/other configs

All configurations are spcified via a config file, whose template is the `config.yaml` file in this repository:

```bash
TODO: contents of config.yaml
```

TODO: edit below

- `<hifi.fastq>`: Input HiFi read dataset.
- `-e`: Threshold of the read depth of contigs. Any contig whose read depth is lower than or equal to this value is filtered out.
- `-T`: Number of threads.

## Output

In each subdirectory with an integer name (e.g. `1/`, `2/`, etc.), 1) input reads are assembled, 2) reads are mapped to the contigs, and 3) reads mapped to low-coverage contigs are filtered out. The filtered reads are used as the input for the next round. The iteration stops when there no longer exist low-coverage contigs.

- `<N-1>/filtered.fastq` (where `N` is the maximum number among the subdirectories): The final filtered read set.
