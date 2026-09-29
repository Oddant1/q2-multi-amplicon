# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model
from qiime2.plugin import SemanticType, ValidationError
from q2_types.sample_data import SampleData


SAMOutput = SemanticType('SAMOutput', variant_of=SampleData.field['type'])

# An alignment record has 11 mandatory fields, optionally followed by tags.
_REQUIRED_FIELD_COUNT = 11

# Zero-based indices of the mandatory fields that must be integers:
# FLAG, POS, MAPQ, PNEXT and TLEN.
_INTEGER_FIELDS = (1, 3, 4, 7, 8)

_FIELD_NAMES = {1: 'FLAG', 3: 'POS', 4: 'MAPQ', 7: 'PNEXT', 8: 'TLEN'}

# SAM bitwise flags for records that do not describe a read's own primary
# alignment: 4 (unmapped), 256 (secondary) and 2048 (supplementary). Skipping
# them means every read is seen exactly once.
_SKIP_FLAGS = 4 | 256 | 2048

# SAM uses '*' where a field is unavailable.
_UNAVAILABLE = '*'


class SAMFormat(model.TextFileFormat):
    """A SAM (Sequence Alignment/Map) file, as written by bowtie2.

    Lines beginning with ``@`` are header lines and are not inspected beyond
    being recognized as headers. Every other line is an alignment record: a
    tab-separated line of at least 11 fields, of which FLAG, POS, MAPQ,
    PNEXT and TLEN must be integers.

    A file containing only headers is valid. bowtie2 writes one when a
    sample has no alignments to report, which happens routinely when it is
    run with ``--no-unal``.

    """

    def _validate_(self, level):
        max_records = {'min': 10, 'max': None}[level]

        with self.open() as fh:
            records_seen = 0
            for line_number, line in enumerate(fh, start=1):
                line = line.rstrip('\n')
                if not line or line.startswith('@'):
                    continue

                self._validate_record(line, line_number)

                records_seen += 1
                if max_records is not None and records_seen >= max_records:
                    break

    def _validate_record(self, line, line_number):
        fields = line.split('\t')

        if len(fields) < _REQUIRED_FIELD_COUNT:
            raise ValidationError(
                'Alignment record on line %d has %d tab-separated field(s), '
                'but SAM requires at least %d.'
                % (line_number, len(fields), _REQUIRED_FIELD_COUNT))

        for index in _INTEGER_FIELDS:
            try:
                int(fields[index])
            except ValueError:
                raise ValidationError(
                    'The %s field of the alignment record on line %d is %r, '
                    'which is not an integer.'
                    % (_FIELD_NAMES[index], line_number, fields[index]))


class SAMDirFmt(model.DirectoryFormat):
    """One SAM file per sample, named after the sample it describes."""

    sams = model.FileCollection(r'.+\.sam', format=SAMFormat)

    @sams.set_path_maker
    def sams_path_maker(self, sample_id):
        return '%s.sam' % sample_id


def _iter_primary_records(sam_path):
    """Yield the fields of each primary alignment record in a SAM file."""
    with open(sam_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('@'):
                continue

            fields = line.split('\t')
            if int(fields[1]) & _SKIP_FLAGS:
                continue

            yield fields
