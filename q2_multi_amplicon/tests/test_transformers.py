# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import gzip
from pathlib import Path

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_types.per_sample_sequences import (
    QIIME1DemuxDirFmt, SingleLanePerSampleSingleEndFastqDirFmt)

from q2_multi_amplicon._types_and_formats import SAMDirFmt

HEADER = '@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:%s\tLN:200\n'


def _record(qname, flag, sequence, quality, reference='ref1'):
    return '\t'.join([qname, str(flag), reference, '1', '42',
                      '%dM' % len(sequence), '*', '0', '0',
                      sequence, quality]) + '\n'


class TransformerTestBase(TestPluginBase):
    package = 'q2_multi_amplicon.tests'

    def _sam_dir(self, records):
        """Write one sample's records to one SAM file per reference."""
        root = Path(self.temp_dir.name) / 'alignments'
        root.mkdir(exist_ok=True)

        by_reference = collections.defaultdict(list)
        for record in records:
            by_reference[record.split('\t')[2]].append(record)
        for reference, reference_records in by_reference.items():
            (root / ('%s.sam' % reference)).write_text(
                HEADER % reference + ''.join(reference_records))
        return str(root)

    def _to_fastq(self, records):
        result = transform(self._sam_dir(records), from_type=SAMDirFmt,
                           to_type=SingleLanePerSampleSingleEndFastqDirFmt)
        observed = {}
        for path in sorted(result.path.glob('*.fastq.gz')):
            with gzip.open(str(path), 'rt') as fh:
                lines = fh.read().splitlines()
            observed[path.name.split('_')[0]] = [
                (lines[i], lines[i + 1], lines[i + 3])
                for i in range(0, len(lines), 4)]
        return result, observed

    def _to_fasta(self, records):
        result = transform(self._sam_dir(records), from_type=SAMDirFmt,
                           to_type=QIIME1DemuxDirFmt)
        lines = (result.path / 'seqs.fna').read_text().splitlines()
        return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines), 2)]


class SAMToSequencesWithQualityTests(TransformerTestBase):
    def test_one_entry_per_reference(self):
        _, observed = self._to_fastq([
            _record('read1', 0, 'ACGTA', 'IIIII'),
            _record('read2', 0, 'TTTTT', 'HHHHH'),
            _record('read3', 0, 'GGCCA', 'JJJJJ', reference='ref2'),
        ])

        self.assertEqual(observed, {
            'ref1': [('@read1', 'ACGTA', 'IIIII'),
                     ('@read2', 'TTTTT', 'HHHHH')],
            'ref2': [('@read3', 'GGCCA', 'JJJJJ')],
        })

    def test_reverse_strand_records_are_restored(self):
        # SAM stores the reverse complement of a read that aligned to the
        # reverse strand, so the read as sequenced is ACGTT/EDCBA.
        _, observed = self._to_fastq([_record('read1', 16, 'AACGT', 'ABCDE')])

        self.assertEqual(observed, {'ref1': [('@read1', 'ACGTT', 'EDCBA')]})

    def test_skips_unmapped_and_non_primary_records(self):
        _, observed = self._to_fastq([
            _record('kept', 0, 'ACGTA', 'IIIII'),
            _record('unmapped', 4, 'TTTTT', 'IIIII'),
            _record('secondary', 256, 'GGGGG', 'IIIII'),
            _record('supplementary', 2048, 'CCCCC', 'IIIII'),
        ])

        self.assertEqual(observed, {'ref1': [('@kept', 'ACGTA', 'IIIII')]})

    def test_writes_manifest_and_metadata(self):
        result, _ = self._to_fastq([
            _record('read1', 0, 'ACGTA', 'IIIII'),
            _record('read2', 0, 'TTTTT', 'IIIII', reference='ref2'),
        ])

        manifest = (result.path / 'MANIFEST').read_text().splitlines()
        self.assertEqual(manifest[0], 'sample-id,filename,direction')
        self.assertEqual([line.split(',')[0] for line in manifest[1:]],
                         ['ref1', 'ref2'])
        self.assertEqual([line.split(',')[2] for line in manifest[1:]],
                         ['forward', 'forward'])
        self.assertIn('phred-offset: 33',
                      (result.path / 'metadata.yml').read_text())

    def test_missing_quality_is_an_error(self):
        with self.assertRaisesRegex(
                ValueError, r"2 read\(s\) aligned to 'ref2' have no "
                            'quality scores'):
            self._to_fastq([
                _record('read1', 0, 'ACGTA', 'IIIII'),
                _record('read2', 0, 'ACGTA', '*', reference='ref2'),
                _record('read3', 0, 'TTTTT', 'IIIII', reference='ref2'),
                _record('read4', 16, 'GGCCA', '*', reference='ref2'),
            ])

    def test_skips_records_without_a_sequence(self):
        _, observed = self._to_fastq([
            _record('read1', 0, 'ACGTA', 'IIIII'),
            _record('read2', 0, '*', '*'),
        ])

        self.assertEqual(observed, {'ref1': [('@read1', 'ACGTA', 'IIIII')]})


class SAMToSequencesTests(TransformerTestBase):
    def test_ids_are_reference_id_and_index(self):
        observed = self._to_fasta([
            _record('read1', 0, 'ACGTA', 'IIIII'),
            _record('read2', 0, 'TTTTT', 'IIIII'),
            _record('read3', 0, 'GGCCA', 'IIIII', reference='ref2'),
        ])

        self.assertEqual(observed, [
            ('ref1_0', 'ACGTA'),
            ('ref1_1', 'TTTTT'),
            ('ref2_0', 'GGCCA'),
        ])

    def test_reverse_strand_records_are_restored(self):
        observed = self._to_fasta([_record('read1', 16, 'AACGT', 'ABCDE')])

        self.assertEqual(observed, [('ref1_0', 'ACGTT')])

    def test_missing_quality_is_not_an_error(self):
        # FASTA carries no quality scores, so a SAM without them transforms
        # cleanly here even though it cannot become SequencesWithQuality.
        observed = self._to_fasta([_record('read1', 0, 'ACGTA', '*')])

        self.assertEqual(observed, [('ref1_0', 'ACGTA')])

    def test_skips_unmapped_and_non_primary_records(self):
        observed = self._to_fasta([
            _record('kept', 0, 'ACGTA', 'IIIII'),
            _record('unmapped', 4, 'TTTTT', 'IIIII'),
            _record('secondary', 256, 'GGGGG', 'IIIII'),
        ])

        self.assertEqual(observed, [('ref1_0', 'ACGTA')])

    def test_skips_records_without_a_sequence(self):
        observed = self._to_fasta([
            _record('read1', 0, '*', '*'),
            _record('read2', 0, 'ACGTA', 'IIIII'),
        ])

        self.assertEqual(observed, [('ref1_0', 'ACGTA')])
