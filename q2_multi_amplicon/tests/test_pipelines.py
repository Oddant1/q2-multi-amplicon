# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip

import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.metadata import ImmutableMetadata
from q2_types.per_sample_sequences import (
    JoinedSequencesWithQuality, SequencesWithQuality,
    SingleLanePerSampleSingleEndFastqDirFmt)
from q2_types.sample_data import SampleData

from q2_multi_amplicon._types_and_formats import SAMDirFmt, SAMOutput

# The primers the test fixtures were built with. See the fixture generation
# notes in tests/data/README.md.
PRIMER_F = 'ACGTACGTACGTACGTACGT'
PRIMER_R = 'TTGGTTGGTTGGTTGGTTGG'

# sample1 holds only ref1 reads, sample3 only ref2 reads, sample2 an even
# split of the two.
EXPECTED_ASSIGNMENTS = {
    'sample1': {'ref1': 10},
    'sample2': {'ref1': 5, 'ref2': 5},
    'sample3': {'ref2': 10},
}


class MultiAmpliconAnalysisTests(TestPluginBase):
    package = 'q2_multi_amplicon.tests'

    def setUp(self):
        super().setUp()
        self.multi_amplicon_analysis = self.plugin.pipelines[
            'multi_amplicon_analysis']
        self.database = qiime2.Artifact.load(
            self.get_data_path('reference-index.qza'))
        self.references = qiime2.Artifact.load(
            self.get_data_path('reference-sequences.qza')).view(pd.Series)
        self.single_end = qiime2.Artifact.load(
            self.get_data_path('single-end-sequences.qza'))
        self.paired_end = qiime2.Artifact.load(
            self.get_data_path('paired-end-sequences.qza'))

    def _assert_expected_classifications(self, classifications):
        # Each sample has a SAM file per reference its reads were assigned
        # to, so the files and their records give the sample's assignments.
        self.assertEqual(list(classifications.keys()),
                         ['sample1', 'sample2', 'sample3'])

        observed = {}
        for sample_id, alignments in classifications.items():
            self.assertEqual(alignments.type, SampleData[SAMOutput])
            # Checks that each file only holds alignments to its reference.
            alignments.validate(level='max')

            observed[sample_id] = {}
            for path in alignments.view(SAMDirFmt).path.glob('*.sam'):
                observed[sample_id][path.stem] = sum(
                    1 for line in path.read_text().splitlines()
                    if not line.startswith('@'))

        self.assertEqual(observed, EXPECTED_ASSIGNMENTS)

    def _assert_expected_table(self, table):
        self.assertEqual(table.type, FeatureTable[Frequency])

        observed = {
            sample_id: {ref: count for ref, count in row.items() if count}
            for sample_id, row in table.view(pd.DataFrame).iterrows()}

        self.assertEqual(observed, EXPECTED_ASSIGNMENTS)

    def _assert_expected_dereplicated(self, results, reference_slice):
        # Every read of a reference is identical once trimmed, so each
        # reference has one feature: the reference sequence itself, in its
        # own orientation, or as much of it as the reads cover.
        references = self.references
        tables = results.dereplicated_tables
        sequences = results.dereplicated_sequences
        self.assertEqual(list(tables.keys()), ['ref1', 'ref2'])
        self.assertEqual(list(sequences.keys()), ['ref1', 'ref2'])

        for reference in 'ref1', 'ref2':
            with self.subTest(reference=reference):
                table = tables[reference]
                self.assertEqual(table.type, FeatureTable[Frequency])
                self.assertEqual(sequences[reference].type,
                                 FeatureData[Sequence])

                seqs = sequences[reference].view(pd.Series)
                self.assertEqual(
                    [str(s) for s in seqs],
                    [str(references[reference])[reference_slice]])

                counts = table.view(pd.DataFrame)[seqs.index[0]]
                self.assertEqual(
                    {sample_id: count for sample_id, count in counts.items()
                     if count},
                    {sample_id: refs[reference]
                     for sample_id, refs in EXPECTED_ASSIGNMENTS.items()
                     if reference in refs})

    def _genotypes(self, results):
        """Map each reference to [sample, read count 1, allele 1, ...] rows."""
        return {
            reference: [[sample_id, *row] for sample_id, row in
                        table.view(qiime2.Metadata).to_dataframe().iterrows()]
            for reference, table in results.genotype_tables.items()}

    def _assert_expected_genotypes(self, results, reference_slice):
        for output, semantic_type in (
                ('alleles', FeatureData[Sequence]),
                ('genotyped_sequences', FeatureData[Sequence]),
                ('genotype_tables', ImmutableMetadata),
                ('allele_frequencies', ImmutableMetadata),
                ('heterozygosity', ImmutableMetadata)):
            collection = getattr(results, output)
            self.assertEqual(list(collection.keys()), ['ref1', 'ref2'])
            for artifact in collection.values():
                self.assertEqual(artifact.type, semantic_type)

        # Each reference has a single allele: the reference itself, or as
        # much of it as the reads cover.
        for reference in 'ref1', 'ref2':
            alleles = results.alleles[reference].view(pd.Series)
            self.assertEqual(
                {allele: str(s) for allele, s in alleles.items()},
                {'A1': str(self.references[reference])[reference_slice]})

        # sample2 has only five reads at each reference, fewer than the
        # default min_coverage of 10.
        self.assertEqual(self._genotypes(results), {
            'ref1': [['sample1', 10, 'A1', 10, 'A1']],
            'ref2': [['sample3', 10, 'A1', 10, 'A1']],
        })

    def _assert_read_lengths(self, sequences, length):
        # Local alignment classifies reads whether or not their primers were
        # trimmed, so the read lengths are what show that trimming happened.
        lengths = set()
        fmt = sequences.view(SingleLanePerSampleSingleEndFastqDirFmt)
        for path in fmt.path.glob('*.fastq.gz'):
            with gzip.open(path, 'rt') as fh:
                lengths.update(
                    len(line.strip()) for i, line in enumerate(fh)
                    if i % 4 == 1)

        self.assertEqual(lengths, {length})

    def _adapter_metadata(self, **columns):
        return qiime2.Metadata(pd.DataFrame(
            {name: [adapter] for name, adapter in columns.items()},
            index=pd.Index(['primer1'], name='id')))

    def test_single_end(self):
        results = self.multi_amplicon_analysis(
            sequences=self.single_end, database=self.database,
            front_f=[PRIMER_F])

        self.assertEqual(results.classified_sequences.type,
                         SampleData[SequencesWithQuality])
        self._assert_expected_classifications(results.classifications)
        self._assert_expected_table(results.table)
        self._assert_expected_dereplicated(results, slice(120))
        self._assert_expected_genotypes(results, slice(120))
        self.assertIn('sample1', results.trim_stats.view(
            qiime2.Metadata).to_dataframe().index)

    def test_paired_end(self):
        # Merging the pairs is what makes the paired-end reads classifiable,
        # so the same counts as the single-end case prove the merge ran.
        results = self.multi_amplicon_analysis(
            sequences=self.paired_end, database=self.database,
            front_f=[PRIMER_F], front_r=[PRIMER_R])

        self.assertEqual(results.classified_sequences.type,
                         SampleData[JoinedSequencesWithQuality])
        self._assert_expected_classifications(results.classifications)
        self._assert_expected_table(results.table)
        self._assert_expected_dereplicated(results, slice(None))
        self._assert_expected_genotypes(results, slice(None))
        self.assertIn('sample1', results.trim_stats.view(
            qiime2.Metadata).to_dataframe().index)

    def test_min_coverage(self):
        results = self.multi_amplicon_analysis(
            sequences=self.single_end, database=self.database,
            front_f=[PRIMER_F], min_coverage=5)

        self.assertEqual(self._genotypes(results), {
            'ref1': [['sample1', 10, 'A1', 10, 'A1'],
                     ['sample2', 5, 'A1', 5, 'A1']],
            'ref2': [['sample2', 5, 'A1', 5, 'A1'],
                     ['sample3', 10, 'A1', 10, 'A1']],
        })

    def test_references_without_enough_coverage_have_no_alleles(self):
        results = self.multi_amplicon_analysis(
            sequences=self.single_end, database=self.database,
            front_f=[PRIMER_F], min_coverage=11)

        self.assertEqual(list(results.dereplicated_tables.keys()),
                         ['ref1', 'ref2'])
        for output in ('alleles', 'genotyped_sequences', 'genotype_tables',
                       'allele_frequencies', 'heterozygosity'):
            self.assertEqual(len(getattr(results, output)), 0)

    def test_reference_alleles_are_matched_per_reference(self):
        ref1 = str(self.references['ref1'])[:120]
        reference_alleles = {'ref1': qiime2.Artifact.import_data(
            'FeatureData[Sequence]',
            pd.Series({'A5': ref1, 'A7': 'ACGTACGTTT'}))}

        results = self.multi_amplicon_analysis(
            sequences=self.single_end, database=self.database,
            front_f=[PRIMER_F], reference_alleles=reference_alleles)

        # ref1 takes its allele names from the reference alleles and keeps
        # all of them, while ref2, which has none, starts from A1.
        self.assertEqual(
            list(results.alleles['ref1'].view(pd.Series).index),
            ['A5', 'A7'])
        self.assertEqual(self._genotypes(results), {
            'ref1': [['sample1', 10, 'A5', 10, 'A5']],
            'ref2': [['sample3', 10, 'A1', 10, 'A1']],
        })

    def test_single_end_metadata(self):
        # The forward-read column is named as it is for paired-end input, so
        # this only trims if the pipeline renamed it for trim_single.
        results = self.multi_amplicon_analysis(
            sequences=self.single_end, database=self.database,
            metadata=self._adapter_metadata(front_f=PRIMER_F))

        self._assert_read_lengths(results.classified_sequences, 120)
        self._assert_expected_classifications(results.classifications)

    def test_single_end_metadata_hyphenated_column(self):
        results = self.multi_amplicon_analysis(
            sequences=self.single_end, database=self.database,
            metadata=self._adapter_metadata(**{'front-f': PRIMER_F}))

        self._assert_read_lengths(results.classified_sequences, 120)

    def test_paired_end_metadata(self):
        results = self.multi_amplicon_analysis(
            sequences=self.paired_end, database=self.database,
            metadata=self._adapter_metadata(
                front_f=PRIMER_F, front_r=PRIMER_R))

        self._assert_read_lengths(results.classified_sequences, 200)
        self._assert_expected_classifications(results.classifications)

    def test_paired_end_metadata_with_params(self):
        results = self.multi_amplicon_analysis(
            sequences=self.paired_end, database=self.database,
            metadata=self._adapter_metadata(front_f=PRIMER_F),
            front_r=[PRIMER_R])

        self._assert_read_lengths(results.classified_sequences, 200)
        self._assert_expected_classifications(results.classifications)

    def test_single_end_rejects_reverse_read_params(self):
        with self.assertRaisesRegex(
                ValueError,
                r'only apply to paired-end sequences.*front_r'):
            self.multi_amplicon_analysis(
                sequences=self.single_end, database=self.database,
                front_f=[PRIMER_F], front_r=[PRIMER_R])

    def test_single_end_rejects_merge_params(self):
        with self.assertRaisesRegex(
                ValueError,
                r'only apply to paired-end sequences.*maxdiffs, minovlen'):
            self.multi_amplicon_analysis(
                sequences=self.single_end, database=self.database,
                front_f=[PRIMER_F], minovlen=25, maxdiffs=0)

    def test_single_end_rejects_reverse_read_metadata(self):
        with self.assertRaisesRegex(
                ValueError,
                r'only apply to paired-end sequences.*'
                r'front-r \(metadata column\)'):
            self.multi_amplicon_analysis(
                sequences=self.single_end, database=self.database,
                metadata=self._adapter_metadata(
                    front_f=PRIMER_F, **{'front-r': PRIMER_R}))
