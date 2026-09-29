# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from pathlib import Path

from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase
from q2_types.sample_data import SampleData

from q2_multi_amplicon._types_and_formats import (
    SAMDirFmt, SAMFormat, SAMOutput)

HEADER = '@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:ref1\tLN:200\n'
RECORD = 'read1\t0\tref1\t1\t42\t20S120M\t*\t0\t0\tACGT\tIIII\n'
REF2_RECORD = RECORD.replace('\tref1\t', '\tref2\t')


class SAMFormatTests(TestPluginBase):
    package = 'q2_multi_amplicon.tests'

    def _format(self, contents):
        path = Path(self.temp_dir.name) / 'alignments.sam'
        path.write_text(contents)
        return SAMFormat(str(path), mode='r')

    def test_valid(self):
        self._format(HEADER + RECORD).validate()

    def test_headers_only_is_valid(self):
        self._format(HEADER).validate()

    def test_optional_tags_are_allowed(self):
        record = RECORD.rstrip('\n') + '\tAS:i:-6\tXN:i:0\n'
        self._format(HEADER + record).validate()

    def test_too_few_fields(self):
        record = 'read1\t0\tref1\t1\t42\n'
        with self.assertRaisesRegex(
                ValidationError, 'line 3 has 5 tab-separated'):
            self._format(HEADER + record).validate()

    def test_non_integer_flag(self):
        record = RECORD.replace('read1\t0\t', 'read1\tforward\t')
        with self.assertRaisesRegex(
                ValidationError, "FLAG field .* line 3 is 'forward'"):
            self._format(HEADER + record).validate()

    def test_non_integer_pos(self):
        record = 'read1\t0\tref1\tfirst\t42\t4M\t*\t0\t0\tACGT\tIIII\n'
        with self.assertRaisesRegex(
                ValidationError, "POS field .* line 3 is 'first'"):
            self._format(HEADER + record).validate()

    def test_min_validation_stops_early(self):
        # The 11th record is malformed, so 'min' validation passes over it
        # while 'max' validation does not.
        contents = HEADER + (RECORD * 10) + 'read11\t0\tref1\n'

        self._format(contents).validate(level='min')

        with self.assertRaisesRegex(ValidationError, 'line 13'):
            self._format(contents).validate(level='max')


class SAMDirFmtTests(TestPluginBase):
    package = 'q2_multi_amplicon.tests'

    def _dirfmt(self, sams):
        root = Path(self.temp_dir.name) / 'alignments'
        root.mkdir()
        for reference, records in sams.items():
            (root / ('%s.sam' % reference)).write_text(
                HEADER + ''.join(records))
        return SAMDirFmt(str(root), mode='r')

    def test_collects_sam_files(self):
        dirfmt = self._dirfmt({'ref1': [RECORD], 'ref2': [REF2_RECORD]})
        dirfmt.validate()

        self.assertEqual(
            sorted(path.stem for path, _ in dirfmt.sams.iter_views(SAMFormat)),
            ['ref1', 'ref2'])

    def test_no_files_is_valid(self):
        # A sample none of whose reads aligned has no alignments to hold.
        self._dirfmt({}).validate()

    def test_record_aligned_to_another_reference(self):
        dirfmt = self._dirfmt({'ref1': [RECORD, REF2_RECORD]})

        with self.assertRaisesRegex(
                ValidationError,
                "record 2 of ref1.sam is aligned to 'ref2'.*only hold "
                "alignments to 'ref1'"):
            dirfmt.validate()

    def test_min_validation_checks_first_records_only(self):
        # The 11th record is aligned elsewhere, so 'min' validation passes
        # over it while 'max' validation does not.
        dirfmt = self._dirfmt({'ref1': [RECORD] * 10 + [REF2_RECORD]})

        dirfmt.validate(level='min')

        with self.assertRaisesRegex(ValidationError, 'record 11 of ref1.sam'):
            dirfmt.validate(level='max')

    def test_path_maker_names_file_after_reference(self):
        dirfmt = SAMDirFmt()

        observed = dirfmt.sams.path_maker(reference='ref1')

        self.assertEqual(observed.name, 'ref1.sam')


class SAMOutputTests(TestPluginBase):
    package = 'q2_multi_amplicon.tests'

    def test_registered_as_sample_data_variant(self):
        self.assertEqual(
            self.plugin.type_fragments['SAMOutput'].fragment, SAMOutput)
        # Raises if SAMOutput is not a valid SampleData field member.
        SampleData[SAMOutput]

    def test_semantic_type_to_format_registration(self):
        self.assertSemanticTypeRegisteredToFormat(
            SampleData[SAMOutput], SAMDirFmt)
