# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import hashlib
import subprocess
import unittest.mock as mock
from pathlib import Path

import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt, QIIME1DemuxDirFmt)

from q2_multi_amplicon._methods import (
    classify_reads, count_classifications, dereplicate_classifications,
    extend_index, _fasta_ids)
from q2_multi_amplicon._types_and_formats import SAMDirFmt


def _rnames(path):
    """The reference that each record of a SAM file is aligned to."""
    return [line.split('\t')[2] for line in path.read_text().splitlines()
            if not line.startswith('@')]


def _references_per_sample(classifications):
    """Map each sample to a count of the references its reads aligned to."""
    observed = {}
    for sample_id, alignments in classifications.items():
        counts = collections.Counter()
        for path in alignments.path.glob('*.sam'):
            counts.update(_rnames(path))
        observed[sample_id] = dict(counts)
    return observed


def _fake_run(observed_cmds, sam=''):
    """Stand in for subprocess.run, recording each command and writing `sam`
    wherever bowtie2 or samtools would have written their output."""
    def fake_run(cmd, **kwargs):
        observed_cmds.append(cmd)
        output_flag = '-S' if cmd[0] == 'bowtie2' else '-o'
        Path(cmd[cmd.index(output_flag) + 1]).write_text(sam)
        return mock.DEFAULT
    return fake_run


class ClassifyReadsTests(TestPluginBase):
    package = 'q2_multi_amplicon.tests'

    def setUp(self):
        super().setUp()
        self.sequences = qiime2.Artifact.load(
            self.get_data_path('single-end-sequences.qza')
        ).view(CasavaOneEightSingleLanePerSampleDirFmt)
        self.database = qiime2.Artifact.load(
            self.get_data_path('reference-index.qza')
        ).view(Bowtie2IndexDirFmt)

    def test_one_result_per_sample(self):
        observed = classify_reads(self.sequences, self.database)

        self.assertEqual(list(observed), ['sample1', 'sample2', 'sample3'])
        for alignments in observed.values():
            self.assertIsInstance(alignments, SAMDirFmt)

    def test_one_sam_per_reference(self):
        observed = classify_reads(self.sequences, self.database)

        self.assertEqual(
            {sample_id: sorted(p.name for p in alignments.path.glob('*.sam'))
             for sample_id, alignments in observed.items()},
            {'sample1': ['ref1.sam'],
             'sample2': ['ref1.sam', 'ref2.sam'],
             'sample3': ['ref2.sam']})

    def test_each_sam_holds_only_its_reference(self):
        observed = classify_reads(self.sequences, self.database)

        for sample_id, alignments in observed.items():
            alignments.validate(level='max')
            for path in alignments.path.glob('*.sam'):
                with self.subTest(sample_id=sample_id, file=path.name):
                    self.assertEqual(set(_rnames(path)), {path.stem})

    def test_assigns_reads_to_expected_references(self):
        # The fixtures were built so that sample1 contains only ref1 reads,
        # sample3 only ref2 reads, and sample2 an even split of the two.
        observed = classify_reads(self.sequences, self.database)

        self.assertEqual(_references_per_sample(observed), {
            'sample1': {'ref1': 10},
            'sample2': {'ref1': 5, 'ref2': 5},
            'sample3': {'ref2': 10},
        })

    def test_headers_describe_only_their_reference(self):
        observed = classify_reads(self.sequences, self.database)

        lines = (observed['sample2'].path / 'ref1.sam').read_text()
        header = [line for line in lines.splitlines() if line.startswith('@')]
        # The @PG lines are bowtie2's and samtools sort's.
        self.assertEqual([line.split('\t')[0] for line in header],
                         ['@HD', '@SQ', '@PG', '@PG'])
        self.assertEqual(header[1], '@SQ\tSN:ref1\tLN:200')

    def test_global_mode(self):
        # bowtie2's end-to-end presets tolerate the untrimmed primer on these
        # reads, so switching modes does not change the assignments.
        observed = classify_reads(self.sequences, self.database, mode='global')

        self.assertEqual(_references_per_sample(observed), {
            'sample1': {'ref1': 10},
            'sample2': {'ref1': 5, 'ref2': 5},
            'sample3': {'ref2': 10},
        })

    def test_builds_expected_command(self):
        observed_cmds = []

        with mock.patch('q2_multi_amplicon._methods.subprocess.run',
                        side_effect=_fake_run(observed_cmds)):
            classify_reads(self.sequences, self.database, n_threads=3,
                           mode='local', sensitivity='very-sensitive')

        # Each sample is aligned, and then its alignments are sorted.
        bowtie2_cmds, sort_cmds = observed_cmds[0::2], observed_cmds[1::2]
        self.assertEqual(len(bowtie2_cmds), 3)
        self.assertEqual(len(sort_cmds), 3)
        for bowtie2_cmd, sort_cmd in zip(bowtie2_cmds, sort_cmds):
            self.assertEqual(bowtie2_cmd[0], 'bowtie2')
            self.assertIn('--very-sensitive-local', bowtie2_cmd)
            self.assertIn('--no-unal', bowtie2_cmd)
            self.assertEqual(bowtie2_cmd[bowtie2_cmd.index('-p') + 1], '3')

            self.assertEqual(sort_cmd[:2], ['samtools', 'sort'])
            self.assertEqual(sort_cmd[-1],
                             bowtie2_cmd[bowtie2_cmd.index('-S') + 1])

        self.assertEqual(
            [Path(cmd[cmd.index('-U') + 1]).name.split('_')[0]
             for cmd in bowtie2_cmds],
            ['sample1', 'sample2', 'sample3'])

    def test_global_mode_command(self):
        observed_cmds = []

        with mock.patch('q2_multi_amplicon._methods.subprocess.run',
                        side_effect=_fake_run(observed_cmds)):
            classify_reads(self.sequences, self.database, mode='global')

        for cmd in observed_cmds[0::2]:
            self.assertIn('--sensitive', cmd)
            self.assertNotIn('--sensitive-global', cmd)

    def test_sample_without_alignments_has_no_files(self):
        # bowtie2 writes a header-only SAM for a sample none of whose reads
        # aligned, which happens routinely under --no-unal.
        with mock.patch('q2_multi_amplicon._methods.subprocess.run',
                        side_effect=_fake_run([], sam=SAM_HEADER)):
            observed = classify_reads(self.sequences, self.database)

        self.assertEqual(
            {sample_id: list(alignments.path.iterdir())
             for sample_id, alignments in observed.items()},
            {'sample1': [], 'sample2': [], 'sample3': []})

    def test_bowtie2_failure_is_reported(self):
        error = subprocess.CalledProcessError(
            returncode=1, cmd=['bowtie2'], stderr='Error: reads file is bad')

        with mock.patch('q2_multi_amplicon._methods.subprocess.run',
                        side_effect=error):
            with self.assertRaisesRegex(
                    RuntimeError,
                    'bowtie2 failed with exit code 1: '
                    'Error: reads file is bad'):
                classify_reads(self.sequences, self.database)


SAM_HEADER = ('@HD\tVN:1.0\tSO:unsorted\n'
              '@SQ\tSN:ref1\tLN:200\n@SQ\tSN:ref2\tLN:200\n')


def _sam_record(qname, flag, reference, sequence='ACGTA'):
    return '\t'.join([qname, str(flag), reference, '1', '42',
                      '%dM' % len(sequence), '*', '0', '0', sequence,
                      'I' * len(sequence)]) + '\n'


class SAMInputTestBase(TestPluginBase):
    package = 'q2_multi_amplicon.tests'

    def _classifications(self, samples):
        """Lay out each sample's records as classify_reads does, with the
        records of each reference in a SAM file of their own."""
        classifications = {}
        for sample_id, records in samples.items():
            root = Path(self.temp_dir.name) / sample_id
            root.mkdir()

            by_reference = collections.defaultdict(list)
            for record in records:
                by_reference[record.split('\t')[2]].append(record)
            for reference, reference_records in by_reference.items():
                (root / ('%s.sam' % reference)).write_text(
                    SAM_HEADER + ''.join(reference_records))

            classifications[sample_id] = SAMDirFmt(str(root), mode='r')
        return classifications

    def _classify_fixtures(self):
        sequences = qiime2.Artifact.load(
            self.get_data_path('single-end-sequences.qza')
        ).view(CasavaOneEightSingleLanePerSampleDirFmt)
        database = qiime2.Artifact.load(
            self.get_data_path('reference-index.qza')
        ).view(Bowtie2IndexDirFmt)
        return classify_reads(sequences, database)


class CountClassificationsTests(SAMInputTestBase):
    def _count(self, alignments):
        return count_classifications(alignments).to_dataframe(
            dense=True).astype(int)

    def test_counts_classify_reads_output(self):
        observed = self._count(self._classify_fixtures())

        pd.testing.assert_frame_equal(observed, pd.DataFrame(
            {'sample1': [10, 0], 'sample2': [5, 5], 'sample3': [0, 10]},
            index=['ref1', 'ref2']))

    def test_counts_only_primary_alignments(self):
        observed = self._count(self._classifications({
            'sample1': [
                _sam_record('read1', 0, 'ref1'),
                _sam_record('read2', 16, 'ref1'),
                _sam_record('unmapped', 4, 'ref1'),
                _sam_record('read1', 256, 'ref2'),
                _sam_record('read2', 2048, 'ref2'),
            ],
        }))

        pd.testing.assert_frame_equal(
            observed, pd.DataFrame({'sample1': [2]}, index=['ref1']))

    def test_sample_without_assignments_is_kept(self):
        observed = self._count(self._classifications({
            'sample1': [_sam_record('read1', 0, 'ref1')],
            'sample2': [],
        }))

        pd.testing.assert_frame_equal(observed, pd.DataFrame(
            {'sample1': [1], 'sample2': [0]}, index=['ref1']))

    def test_no_assignments_in_any_sample(self):
        observed = count_classifications(self._classifications({
            'sample1': [],
            'sample2': [_sam_record('unmapped', 4, 'ref1')],
        }))

        self.assertEqual(list(observed.ids()), ['sample1', 'sample2'])
        self.assertEqual(list(observed.ids(axis='observation')), [])


def _md5(sequence):
    return hashlib.md5(sequence.encode()).hexdigest()


class DereplicateClassificationsTests(SAMInputTestBase):
    def _dereplicate(self, alignments):
        """Map each reference to its table and its sequences as strings."""
        tables, sequences = dereplicate_classifications(alignments)
        self.assertEqual(list(tables), list(sequences))
        return {reference: (
                    tables[reference].to_dataframe(dense=True).astype(int),
                    sequences[reference].astype(str))
                for reference in tables}

    def test_dereplicates_classify_reads_output(self):
        observed = self._dereplicate(self._classify_fixtures())

        self.assertEqual(list(observed), ['ref1', 'ref2'])
        # Every read of a reference is identical in the fixtures, so each
        # reference has a single feature holding all of its reads.
        for reference, counts in (('ref1', [10, 5, 0]),
                                  ('ref2', [0, 5, 10])):
            table, sequences = observed[reference]
            feature_id, = table.index
            pd.testing.assert_frame_equal(table, pd.DataFrame(
                [counts], index=[feature_id],
                columns=['sample1', 'sample2', 'sample3']))
            self.assertEqual(list(sequences.index), [feature_id])
            self.assertEqual(feature_id, _md5(sequences[feature_id]))
            # The untrimmed 20 bp primer plus 120 bp of the reference.
            self.assertEqual(len(sequences[feature_id]), 140)

    def test_distinct_sequences_are_separate_features(self):
        observed = self._dereplicate(self._classifications({
            'sample1': [
                _sam_record('read1', 0, 'ref1', 'AAAAA'),
                _sam_record('read2', 0, 'ref1', 'CCCCC'),
                _sam_record('read3', 0, 'ref1', 'CCCCC'),
                _sam_record('read4', 0, 'ref1', 'TTTTT'),
            ],
            'sample2': [
                _sam_record('read5', 0, 'ref1', 'AAAAA'),
                _sam_record('read6', 0, 'ref1', 'AAAAA'),
                _sam_record('read7', 0, 'ref1', 'GGGGG'),
            ],
        }))

        table, sequences = observed['ref1']
        # Most abundant first, and GGGGG before TTTTT because they tie.
        expected_order = ['AAAAA', 'CCCCC', 'GGGGG', 'TTTTT']
        self.assertEqual(list(sequences), expected_order)
        self.assertEqual(list(sequences.index),
                         [_md5(s) for s in expected_order])
        pd.testing.assert_frame_equal(table, pd.DataFrame(
            {'sample1': [1, 2, 0, 1], 'sample2': [2, 0, 1, 0]},
            index=sequences.index))

    def test_reverse_strand_reads_collapse_with_forward_reads(self):
        # SAM already stores a reverse-strand read as its reverse
        # complement, so both records here are the same sequence.
        observed = self._dereplicate(self._classifications({
            'sample1': [_sam_record('read1', 0, 'ref1', 'ACGTT'),
                        _sam_record('read2', 16, 'ref1', 'ACGTT')],
        }))

        table, sequences = observed['ref1']
        self.assertEqual(list(sequences), ['ACGTT'])
        self.assertEqual(table.values.tolist(), [[2]])

    def test_one_table_per_reference_with_every_sample(self):
        observed = self._dereplicate(self._classifications({
            'sample1': [_sam_record('read1', 0, 'ref1', 'AAAAA')],
            'sample2': [_sam_record('read2', 0, 'ref2', 'CCCCC')],
        }))

        self.assertEqual(list(observed), ['ref1', 'ref2'])
        self.assertEqual(observed['ref1'][0].values.tolist(), [[1, 0]])
        self.assertEqual(observed['ref2'][0].values.tolist(), [[0, 1]])
        for table, _ in observed.values():
            self.assertEqual(list(table.columns), ['sample1', 'sample2'])

    def test_uses_only_primary_alignments(self):
        observed = self._dereplicate(self._classifications({
            'sample1': [
                _sam_record('read1', 0, 'ref1', 'AAAAA'),
                _sam_record('unmapped', 4, 'ref1', 'CCCCC'),
                _sam_record('read1', 256, 'ref2', 'AAAAA'),
                _sam_record('read1', 2048, 'ref1', 'GGGGG'),
            ],
        }))

        self.assertEqual(list(observed), ['ref1'])
        self.assertEqual(list(observed['ref1'][1]), ['AAAAA'])

    def test_no_assignments_in_any_sample(self):
        tables, sequences = dereplicate_classifications(self._classifications({
            'sample1': [],
            'sample2': [_sam_record('unmapped', 4, 'ref1')],
        }))

        self.assertEqual((tables, sequences), ({}, {}))


def _index_entries(database):
    base = str(database.path / database.get_basename())
    names = subprocess.run(['bowtie2-inspect', '-n', base],
                           check=True, capture_output=True, text=True)
    return names.stdout.split()


class ExtendIndexTests(TestPluginBase):
    package = 'q2_multi_amplicon.tests'

    def setUp(self):
        super().setUp()
        self.database_artifact = qiime2.Artifact.load(
            self.get_data_path('reference-index.qza'))
        self.database = self.database_artifact.view(Bowtie2IndexDirFmt)
        self.sequences = qiime2.Artifact.load(
            self.get_data_path('single-end-sequences.qza')
        ).view(QIIME1DemuxDirFmt)

    def test_combines_existing_and_new_sequences(self):
        observed = extend_index(self.sequences, self.database)

        entries = _index_entries(observed)
        # The two references already in the index, plus the 30 reads.
        self.assertEqual(len(entries), 32)
        self.assertEqual(entries[:2], ['ref1', 'ref2'])
        self.assertIn('sample1_0', entries)
        # q2-types numbers the sequences continuously across the file rather
        # than restarting per sample, so sample3's last read is _29.
        self.assertIn('sample3_29', entries)

    def test_original_index_is_untouched(self):
        extend_index(self.sequences, self.database)

        self.assertEqual(_index_entries(self.database), ['ref1', 'ref2'])

    def test_accepts_both_input_types(self):
        action = self.plugin.methods['extend_index']

        for fixture in ('single-end-sequences.qza',
                        'single-end-sequences-no-quality.qza'):
            with self.subTest(fixture=fixture):
                sequences = qiime2.Artifact.load(
                    self.get_data_path(fixture))
                observed, = action(sequences=sequences,
                                   database=self.database_artifact)

                self.assertEqual(
                    len(_index_entries(observed.view(Bowtie2IndexDirFmt))), 32)

    def test_rejects_duplicate_ids(self):
        def fake_run(cmd, **kwargs):
            if cmd[0] == 'bowtie2-inspect':
                kwargs['stdout'].write('>sample1_0\nACGTACGT\n'
                                       '>sample2_10\nTTTTTTTT\n')
            return mock.DEFAULT

        with mock.patch('q2_multi_amplicon._methods.subprocess.run',
                        side_effect=fake_run):
            with self.assertRaisesRegex(
                    ValueError,
                    'already in the index.*sample1_0, sample2_10'):
                extend_index(self.sequences, self.database)

    def test_builds_new_index_from_combined_fasta(self):
        combined = {}

        def fake_run(cmd, **kwargs):
            if cmd[0] == 'bowtie2-inspect':
                kwargs['stdout'].write('>ref1\nACGTACGT\n')
            else:
                combined['cmd'] = cmd
                combined['fasta'] = open(cmd[-2]).read()
            return mock.DEFAULT

        with mock.patch('q2_multi_amplicon._methods.subprocess.run',
                        side_effect=fake_run):
            extend_index(self.sequences, self.database, n_threads=3)

        self.assertEqual(combined['cmd'][0], 'bowtie2-build')
        self.assertEqual(combined['cmd'][combined['cmd'].index('--threads')
                                         + 1], '3')
        self.assertTrue(combined['cmd'][-1].endswith('/db'))

        headers = [line for line in combined['fasta'].splitlines()
                   if line.startswith('>')]
        self.assertEqual(headers[0], '>ref1')
        self.assertEqual(headers[1], '>sample1_0')
        self.assertEqual(len(headers), 31)

    def test_bowtie2_inspect_failure_is_reported(self):
        error = subprocess.CalledProcessError(
            returncode=1, cmd=['bowtie2-inspect'],
            stderr='Error: could not open index')

        with mock.patch('q2_multi_amplicon._methods.subprocess.run',
                        side_effect=error):
            with self.assertRaisesRegex(
                    RuntimeError,
                    'bowtie2-inspect failed with exit code 1: '
                    'Error: could not open index'):
                extend_index(self.sequences, self.database)


class FastaIdTests(TestPluginBase):
    package = 'q2_multi_amplicon.tests'

    def _fasta(self, contents):
        path = Path(self.temp_dir.name) / 'seqs.fasta'
        path.write_text(contents)
        return str(path)

    def test_ids_stop_at_whitespace(self):
        observed = _fasta_ids(self._fasta(
            '>ref1 the first reference\nACGT\n>ref2\tsecond\nTTTT\n'))

        self.assertEqual(observed, {'ref1', 'ref2'})

    def test_empty_header_is_skipped(self):
        observed = _fasta_ids(self._fasta('>\nACGT\n>ref1\nTTTT\n'))

        self.assertEqual(observed, {'ref1'})
