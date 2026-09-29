# q2-multi-amplicon

A [QIIME 2](https://qiime2.org) plugin [developed](https://develop.qiime2.org) by Oddant1 (oddant1@hotmail.com). 🔌

## Installation instructions

**The following instructions are intended to be a starting point** and should be replaced when `q2-multi-amplicon` is ready to share with others.
They will enable you to install the most recent *development* version of `q2-multi-amplicon`.
Remember that *release* versions should be used for all "real" work (i.e., where you're not testing or prototyping) - if there aren't instructions for installing a release version of this plugin, it is probably not yet intended for use in practice.

### Install Prerequisites

[Miniconda](https://conda.io/miniconda.html) provides the `conda` environment and package manager, and is currently the only supported way to install QIIME 2.
Follow the instructions for downloading and installing Miniconda.

After installing Miniconda and opening a new terminal, make sure you're running the latest version of `conda`:

```bash
conda update conda
```

###  Install development version of `q2-multi-amplicon`

Next, you need to get into the top-level `q2-multi-amplicon` directory.
If you already have this (e.g., because you just created the plugin), this may be as simple as running `cd q2-multi-amplicon`.
If not, you'll need the `q2-multi-amplicon` directory on your computer.
How you do that will differ based on how the package is shared, and ideally the developer will update these instructions to be more specific (remember, these instructions are intended to be a starting point).
For example, if it's maintained in a GitHub repository, you can achieve this by [cloning the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).
Once you have the directory on your computer, change (`cd`) into it.

If you're in a conda environment, deactivate it by running `conda deactivate`.


Then, follow the install instructions below, based on your machine's architecture:

<details>
<summary><strong>🍏&nbsp;Apple Silicon (ARM)</strong></summary>
<p>&nbsp;</p>

Start by creating a new conda environment:

```shell
CONDA_SUBDIR=osx-64 conda env create -n q2-multi-amplicon-dev --file ./environment-files/q2-multi-amplicon-qiime2-qiime2-dev.yml
```

After this completes, activate the new environment you created by running:

```shell
conda activate q2-multi-amplicon-dev
```

Once this new environment has been activated, update your conda config to set the subdir to osx-64:

```shell
conda config --env --set subdir osx-64
```

Finally, run:

```shell
make install
```
</details>

<details>
<summary><strong>🛠&nbsp;All other architectures (Apple Intel, Linux, WSL)</strong></summary>
<p>&nbsp;</p>

Start by creating a new conda environment:

```shell
conda env create -n q2-multi-amplicon-dev --file ./environment-files/q2-multi-amplicon-qiime2-qiime2-dev.yml
```

After this completes, activate the new environment you created by running:

```shell
conda activate q2-multi-amplicon-dev
```

Finally, run:

```shell
make install
```
</details>

## Testing and using the most recent development version of `q2-multi-amplicon`

After completing the install steps above, confirm that everything is working as expected by running:

```shell
make test
```

You should get a report that tests were run, and you should see that all tests passed and none failed.
It's usually ok if some warnings are reported.

If all of the tests pass, you're ready to use the plugin.
Start by making QIIME 2's command line interface aware of `q2-multi-amplicon` by running:

```shell
qiime dev refresh-cache
```

You should then see the plugin in the list of available plugins if you run:

```shell
qiime info
```

You should be able to review the help text by running:

```shell
qiime multi-amplicon --help
```

Have fun! 😎

## Usage

### `multi-amplicon-analysis`

Trims primers and adapters from demultiplexed reads with
[cutadapt](https://cutadapt.readthedocs.io), merges read pairs with
[vsearch](https://github.com/torognes/vsearch) when the input is paired-end,
assigns the resulting reads, per sample, to the reference sequences in a
[bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) index, counts the
reads assigned to each reference in each sample — in total, and per distinct
read sequence — and finally calls the alleles and genotypes of every sample at
each reference with
[q2-real-allele](#calling-alleles-with-q2-real-allele)'s `call-alleles`.
q2-real-allele must be installed alongside this plugin.

First build an index from your reference amplicon sequences:

```shell
qiime quality-control bowtie2-build \
  --i-sequences reference-sequences.qza \
  --o-database reference-index.qza
```

Then run the pipeline on single-end reads:

```shell
qiime multi-amplicon multi-amplicon-analysis \
  --i-sequences single-end-sequences.qza \
  --i-database reference-index.qza \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --output-dir results
```

...or on paired-end reads, which adds the merging step:

```shell
qiime multi-amplicon multi-amplicon-analysis \
  --i-sequences paired-end-sequences.qza \
  --i-database reference-index.qza \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r GGACTACNVGGGTWTCTAAT \
  --output-dir results
```

The adapters can also be given as a metadata file in place of the adapter
parameters, which is handy when trimming the primers of several amplicons at
once. Each column is named after one of the adapter parameters (`adapter_f`,
`front_f`, `anywhere_f`, `adapter_r`, `front_r`, `anywhere_r`) and holds one
adapter per row:

```
id	front_f	front_r
16S-V4	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT
ITS1	CTTGGTCATTTAGAGGAAGTAA	GCTGCGTTCTTCATCGATGC
```

```shell
qiime multi-amplicon multi-amplicon-analysis \
  --i-sequences paired-end-sequences.qza \
  --i-database reference-index.qza \
  --m-metadata-file primers.tsv \
  --output-dir results
```

Metadata columns and adapter parameters can be mixed, as long as the same
adapter isn't given both ways (for example, a `front_f` column and
`--p-front-f`).

The pipeline writes these outputs to `results/`:

- `classified_sequences.qza`: the reads that went into the classification —
  trimmed, and merged if the input was paired-end. Pairs that vsearch could not
  merge are not classified and are not returned.
- `classifications.qza`: a `SampleData[SAMOutput]` artifact with one SAM file
  per sample, named after that sample, recording which reference sequence each
  read was assigned to.
- `table.qza`: a `FeatureTable[Frequency]` counting those assignments. Its
  features are the IDs of the reference sequences, and each value is the
  number of a sample's reads that were assigned to that reference (see
  [`count-classifications`](#count-classifications)).
- `dereplicated_tables/` and `dereplicated_sequences/`: collections with one
  artifact per reference sequence, holding the distinct read sequences
  assigned to that reference and how many reads had each one in each sample
  (see [`dereplicate-classifications`](#dereplicate-classifications)).
- `alleles/`, `genotyped_sequences/`, `genotype_tables/`,
  `allele_frequencies/` and `heterozygosity/`: the outputs of
  `qiime real-allele call-alleles` for each reference sequence, again one
  artifact per reference (see
  [Calling alleles](#calling-alleles-with-q2-real-allele)).
- `trim_stats.qza`: the per-sample summary of the cutadapt trimming step.

The reverse-read parameters (`--p-adapter-r`, `--p-front-r`,
`--p-anywhere-r`, `--p-reverse-cut`), their metadata columns, and the merging
parameters (`--p-minovlen`, `--p-maxdiffs`, `--p-maxee`, and so on) only
apply to paired-end input; passing any of them alongside single-end reads is
an error rather than a silent no-op.

### `classify-reads`

The classification step is also registered on its own, so it can be run
against single-end or already-merged reads without the trimming and merging
steps:

```shell
qiime multi-amplicon classify-reads \
  --i-sequences merged-sequences.qza \
  --i-database reference-index.qza \
  --o-classifications classifications.qza
```

Reads that did not align to any reference are left out of the output
entirely.

### `count-classifications`

Turns the per-sample alignments into a feature table, so the assignments can
be used with the rest of QIIME 2 (`qiime feature-table summarize`, diversity
analyses, and so on):

```shell
qiime multi-amplicon count-classifications \
  --i-classifications classifications.qza \
  --o-table table.qza
```

The features of the table are the IDs of the reference sequences in the index,
and each value is the number of a sample's reads that were assigned to that
reference. Only each read's primary alignment is counted, so a read is counted
at most once. Every sample gets a column, including one in which no read was
assigned to any reference; such a sample has a count of zero throughout.
References that no read was assigned to do not appear in the table.

### `dereplicate-classifications`

Splits the alignments up by the reference sequence each read was assigned to,
and collapses identical reads within each reference into one feature:

```shell
qiime multi-amplicon dereplicate-classifications \
  --i-classifications classifications.qza \
  --o-dereplicated-tables dereplicated-tables \
  --o-dereplicated-sequences dereplicated-sequences
```

Both outputs are collections keyed by reference ID, saved as a directory with
one artifact per reference (for example `dereplicated-tables/ref1.qza` and
`dereplicated-sequences/ref1.qza`). For each reference:

- the `FeatureData[Sequence]` holds every distinct read sequence assigned to
  it, most abundant first, each identified by the MD5 hash of its sequence;
- the `FeatureTable[Frequency]` gives the number of reads with each of those
  sequences in each sample. Every sample gets a column, including one with no
  reads assigned to that reference.

SAM stores a read that aligned to the reverse strand as its reverse
complement, so a reference's sequences all share its orientation and a
sequence is counted together whichever strand its reads came from. Only each
read's primary alignment is used, and references that no read was assigned to
are left out.

### Calling alleles with `q2-real-allele`

`multi-amplicon-analysis` ends by running `qiime real-allele call-alleles`
once for each reference sequence, since alleles are named per locus, on that
reference's dereplicated table and sequences. Each allele output is a
collection keyed by reference ID, so `results/genotype_tables/ref1.qza` holds
the genotypes at `ref1`, and `results/alleles/ref1.qza` its allele sequences.

`--p-min-coverage` (default 10) and `--p-proportion` (default 0.5) are passed
on to `call-alleles`. A reference sequence at which no sample has a distinct
sequence with `--p-min-coverage` reads cannot be genotyped, so it is left out
of the allele outputs rather than failing the run; its reads are still in the
other outputs.

To name alleles consistently across runs, pass an earlier run's `alleles`
collection as `--i-reference-alleles`. Each reference sequence is matched with
the entry of the same name, and one without an entry gets new allele names
starting from A1:

```shell
qiime multi-amplicon multi-amplicon-analysis \
  --i-sequences paired-end-sequences-batch-2.qza \
  --i-database reference-index.qza \
  --i-reference-alleles results/alleles \
  --m-metadata-file primers.tsv \
  --output-dir results-batch-2
```

Because the pipeline has already trimmed the primers and put every read in the
orientation of its reference, there is no need for the manual orientation and
primer-removal steps of the standalone `Real_Allele_Pipeline.py`.

To call alleles yourself — for instance, from the output of
[`dereplicate-classifications`](#dereplicate-classifications) — run
`call-alleles` once per reference, passing the matching table and sequences:

```shell
for table in dereplicated-tables/*.qza; do
  locus=$(basename "$table" .qza)
  qiime real-allele call-alleles \
    --i-sequences "dereplicated-sequences/$locus.qza" \
    --i-table "$table" \
    --output-dir "alleles-$locus"
done
```

### `extend-index`

Adds sequences to an existing Bowtie2 index — for instance, folding a new
amplicon region's references into the index you already classify against:

```shell
qiime multi-amplicon extend-index \
  --i-sequences new-region-sequences.qza \
  --i-database reference-index.qza \
  --o-extended-database extended-index.qza
```

`--i-sequences` accepts either `SampleData[SequencesWithQuality]` or
`SampleData[Sequences]`. The existing index's references are written back out
as FASTA with `bowtie2-inspect`, the new sequences are appended, and a fresh
index is built from the combination — the index you passed in is not
modified. If any incoming sequence has an ID that the index already uses, the
action fails rather than building an index in which that name is ambiguous.

### Working with `SampleData[SAMOutput]`

The SAM alignments carry the classified reads themselves, so a
`SampleData[SAMOutput]` artifact can be viewed as either of the per-sample
sequence types:

```python
import qiime2
from q2_types.per_sample_sequences import (
    QIIME1DemuxDirFmt, SingleLanePerSampleSingleEndFastqDirFmt)

classifications = qiime2.Artifact.load('classifications.qza')

# SampleData[SequencesWithQuality] — per-sample FASTQ
classifications.view(SingleLanePerSampleSingleEndFastqDirFmt)

# SampleData[Sequences] — a single seqs.fna with '<sample-id>_<n>' headers
classifications.view(QIIME1DemuxDirFmt)
```

Both transformers use `samtools view` to keep only each read's own primary
alignment, so a read appears at most once, and to undo the reverse complement
that SAM stores for reads that aligned to the reverse strand. Viewing as
`SampleData[SequencesWithQuality]` requires the SAM to carry quality scores;
bowtie2 preserves them from the input FASTQ, but a SAM whose `QUAL` field is
`*` can only be viewed as `SampleData[Sequences]`.

## About

The `q2-multi-amplicon` Python package was [created from a template](https://develop.qiime2.org/en/stable/plugins/tutorials/create-from-template.html).
To learn more about `q2-multi-amplicon`, refer to the [project website](https://example.com/q2-multi-amplicon).
To learn how to use QIIME 2, refer to the [QIIME 2 User Documentation](https://docs.qiime2.org).
To learn QIIME 2 plugin development, refer to [*Developing with QIIME 2*](https://develop.qiime2.org).

`q2-multi-amplicon` is a QIIME 2 community plugin, meaning that it is not necessarily developed and maintained by the developers of QIIME 2.
Please be aware that because community plugins are developed by the QIIME 2 developer community, and not necessarily the QIIME 2 developers themselves, some may not be actively maintained or compatible with current release versions of the QIIME 2 distributions.
More information on development and support for community plugins can be found [here](https://library.qiime2.org).
If you need help with a community plugin, first refer to the [project website](https://example.com/q2-multi-amplicon).
If that page doesn't provide information on how to get help, or you need additional help, head to the [Community Plugins category](https://forum.qiime2.org/c/community-contributions/community-plugins/14) on the QIIME 2 Forum where the QIIME 2 developers will do their best to help you.
