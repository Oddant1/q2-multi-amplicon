# Test data

These fixtures back the `classify_reads`, `count_classifications` and
`multi_amplicon_analysis` tests. They describe two 200 bp reference
"amplicons" (`ref1` and `ref2`) and three samples built from them:

| sample  | ref1 reads | ref2 reads |
| ------- | ---------- | ---------- |
| sample1 | 10         | 0          |
| sample2 | 5          | 5          |
| sample3 | 0          | 10         |

- `reference-sequences.qza` — `FeatureData[Sequence]`, the two references.
- `reference-index.qza` — `Bowtie2Index`, built from those references with
  `qiime quality-control bowtie2-build`.
- `single-end-sequences.qza` — `SampleData[SequencesWithQuality]`. Each read
  is the 20 bp forward primer `ACGTACGTACGTACGTACGT` followed by the first
  120 bp of its reference.
- `single-end-sequences-no-quality.qza` — `SampleData[Sequences]`, the same
  reads as above without quality scores, used to check that `extend_index`
  accepts both per-sample sequence types. Produced by viewing
  `single-end-sequences.qza` as `QIIME1DemuxDirFmt` and re-importing it.
- `paired-end-sequences.qza` — `SampleData[PairedEndSequencesWithQuality]`.
  Forward reads are as above; reverse reads are the 20 bp reverse primer
  `TTGGTTGGTTGGTTGGTTGG` followed by the reverse complement of the last
  120 bp of the reference, leaving a 40 bp overlap so that every pair
  merges into the full 200 bp reference.

All quality scores are `I` (Q40), so nothing is lost to quality filtering
and the expected counts are exact.

## Regenerating

Run the following in an environment with `q2-multi-amplicon` installed, from
this directory.

```python
import gzip
import random
from pathlib import Path

import qiime2
from qiime2.plugins import quality_control

PRIMER_F = 'ACGTACGTACGTACGTACGT'
PRIMER_R = 'TTGGTTGGTTGGTTGGTTGG'
REF_LEN, READ_LEN, REV_START = 200, 120, 80
COMPOSITION = {
    'sample1': {'ref1': 10, 'ref2': 0},
    'sample2': {'ref1': 5, 'ref2': 5},
    'sample3': {'ref1': 0, 'ref2': 10},
}

rng = random.Random(42)
refs = {name: ''.join(rng.choice('ACGT') for _ in range(REF_LEN))
        for name in ('ref1', 'ref2')}


def revcomp(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def fastq_record(rid, seq):
    return '@%s\n%s\n+\n%s\n' % (rid, seq, 'I' * len(seq))


def write_reads(directory, paired):
    directory.mkdir(parents=True, exist_ok=True)
    for i, (sample, counts) in enumerate(COMPOSITION.items()):
        fwd, rev = [], []
        for ref_name, n in counts.items():
            ref = refs[ref_name]
            for j in range(n):
                rid = '%s_%s_%d' % (sample, ref_name, j)
                fwd.append(fastq_record(rid, PRIMER_F + ref[:READ_LEN]))
                if paired:
                    rev.append(fastq_record(
                        rid, PRIMER_R + revcomp(ref[REV_START:])))
        stem = '%s_S%d_L001' % (sample, i + 1)
        with gzip.open(directory / ('%s_R1_001.fastq.gz' % stem), 'wt') as fh:
            fh.write(''.join(fwd))
        if paired:
            with gzip.open(directory / ('%s_R2_001.fastq.gz' % stem),
                           'wt') as fh:
                fh.write(''.join(rev))


fasta = Path('refs.fasta')
fasta.write_text(''.join('>%s\n%s\n' % kv for kv in refs.items()))
ref_seqs = qiime2.Artifact.import_data('FeatureData[Sequence]', str(fasta))
ref_seqs.save('reference-sequences.qza')
index, = quality_control.actions.bowtie2_build(sequences=ref_seqs)
index.save('reference-index.qza')

write_reads(Path('se'), paired=False)
qiime2.Artifact.import_data(
    'SampleData[SequencesWithQuality]', 'se',
    view_type='CasavaOneEightSingleLanePerSampleDirFmt'
).save('single-end-sequences.qza')

write_reads(Path('pe'), paired=True)
qiime2.Artifact.import_data(
    'SampleData[PairedEndSequencesWithQuality]', 'pe',
    view_type='CasavaOneEightSingleLanePerSampleDirFmt'
).save('paired-end-sequences.qza')

from q2_types.per_sample_sequences import QIIME1DemuxDirFmt

single_end = qiime2.Artifact.load('single-end-sequences.qza')
demux = single_end.view(QIIME1DemuxDirFmt)
qiime2.Artifact.import_data(
    'SampleData[Sequences]', str(demux.path)
).save('single-end-sequences-no-quality.qza')
```
