# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip

import pandas as pd
import yaml

from q2_types.per_sample_sequences import (
    FastqManifestFormat, QIIME1DemuxDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt, YamlFormat)

from q2_multi_amplicon.plugin_setup import plugin
from q2_multi_amplicon._types_and_formats import (
    SAMDirFmt, SAMFormat, _UNAVAILABLE, _iter_primary_records)

# SAM stores the reverse complement of the original read when the read
# aligned to the reverse strand (flag 16), so that has to be undone to
# recover the sequence as it was read.
_REVERSE_STRAND = 16

_COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVNacgtryswkmbdhvn',
                            'TGCAYRSWMKVHDBNtgcayrswmkvhdbn')


@plugin.register_transformer
def _1(dirfmt: SAMDirFmt) -> SingleLanePerSampleSingleEndFastqDirFmt:
    result = SingleLanePerSampleSingleEndFastqDirFmt()
    manifest_data = []

    for index, (sample_id, records) in enumerate(_iter_samples(dirfmt)):
        output_fp = result.sequences.path_maker(
            sample_id=sample_id,
            # These are not used internally by QIIME 2, so their values
            # don't matter beyond keeping the filenames unique.
            barcode_id=index, lane_number=1, read_number=1)

        with gzip.open(str(output_fp), 'wt') as fh:
            for read_id, sequence, quality in records:
                if quality is None:
                    raise ValueError(
                        'The alignment record for read %r in sample %r has '
                        'no quality scores, so it cannot be represented as '
                        'SampleData[SequencesWithQuality]. View this data as '
                        'SampleData[Sequences] instead.'
                        % (read_id, sample_id))
                fh.write('@%s\n%s\n+\n%s\n' % (read_id, sequence, quality))

        manifest_data.append([sample_id, output_fp.name, 'forward'])

    manifest = FastqManifestFormat()
    pd.DataFrame(
        manifest_data, columns=manifest.EXPECTED_HEADER
    ).to_csv(str(manifest), index=False)
    result.manifest.write_data(manifest, FastqManifestFormat)

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': 33}))
    result.metadata.write_data(metadata, YamlFormat)

    return result


@plugin.register_transformer
def _2(dirfmt: SAMDirFmt) -> QIIME1DemuxDirFmt:
    result = QIIME1DemuxDirFmt()

    with open(str(result.path / 'seqs.fna'), 'w') as fh:
        for sample_id, records in _iter_samples(dirfmt):
            # QIIME 1 demultiplexed FASTA identifies a sequence by
            # '<sample-id>_<seq-id>', where seq-id is unique within the
            # sample. The read's own name is not part of the format.
            for seq_id, (_, sequence, _) in enumerate(records):
                fh.write('>%s_%d\n%s\n' % (sample_id, seq_id, sequence))

    return result


def _iter_samples(dirfmt):
    for path, sam in sorted(dirfmt.sams.iter_views(SAMFormat)):
        yield path.stem, list(_iter_primary_alignments(str(sam)))


def _iter_primary_alignments(sam_path):
    for fields in _iter_primary_records(sam_path):
        read_id, sequence, quality = fields[0], fields[9], fields[10]
        if sequence == _UNAVAILABLE:
            continue
        if quality == _UNAVAILABLE:
            quality = None

        if int(fields[1]) & _REVERSE_STRAND:
            sequence = sequence.translate(_COMPLEMENT)[::-1]
            if quality is not None:
                quality = quality[::-1]

        yield read_id, sequence, quality
