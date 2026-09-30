# q2-multi-amplicon

## Installation instructions

```shell
conda env create -n q2-multi-amplicon-dev -f ./environment-files/q2-multi-amplicon-qiime2-tiny-2026.7.yml
conda activate q2-multi-amplicon-dev
```

 <!-- I asked Claude to wtite a README and it did this... not sure how I feel about it's very long, but it seems potentially useful -->

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
- `classifications/`: a collection with one `SampleData[SAMOutput]` artifact
  per sample (for example `classifications/sample1.qza`), recording which
  reference sequence each of that sample's reads was assigned to. Within each
  sample's artifact, the alignments to each reference are a separate SAM file
  named after that reference (see [`classify-reads`](#classify-reads)).
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
  --o-classifications classifications
```

The output is a collection keyed by sample ID, saved as a directory with one
`SampleData[SAMOutput]` artifact per sample. Each sample's artifact holds one
SAM file per reference sequence that its reads were assigned to, named after
that reference, so `classifications/sample1.qza` might hold `ref1.sam` and
`ref2.sam`. Every record in a file is an alignment to the reference the file
is named after, and each file's header lists only that reference.

Reads that did not align to any reference are left out of the output
entirely. A sample none of whose reads aligned is still in the collection, as
an artifact with no SAM files.

### `count-classifications`

Turns the per-sample alignments into a feature table, so the assignments can
be used with the rest of QIIME 2 (`qiime feature-table summarize`, diversity
analyses, and so on):

```shell
qiime multi-amplicon count-classifications \
  --i-classifications classifications \
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
  --i-classifications classifications \
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

Each `SampleData[SAMOutput]` artifact holds the alignments of one sample, one
SAM file per reference. The SAM alignments carry the classified reads
themselves, so a sample's artifact can be viewed as either of the per-sample
sequence types, with one entry per reference — its reads grouped by the
reference they were assigned to:

```python
import qiime2
from q2_types.per_sample_sequences import (
    QIIME1DemuxDirFmt, SingleLanePerSampleSingleEndFastqDirFmt)

classifications = qiime2.ResultCollection.load('classifications')
sample1 = classifications['sample1']

# SampleData[SequencesWithQuality] — one FASTQ per reference
sample1.view(SingleLanePerSampleSingleEndFastqDirFmt)

# SampleData[Sequences] — a single seqs.fna with '<reference-id>_<n>' headers
sample1.view(QIIME1DemuxDirFmt)
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
