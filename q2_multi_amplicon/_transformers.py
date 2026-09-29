# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import yaml

from q2_types.per_sample_sequences import (
    FastqManifestFormat, QIIME1DemuxDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt, YamlFormat)

from q2_multi_amplicon.plugin_setup import plugin
from q2_multi_amplicon._methods import _run_command
from q2_multi_amplicon._types_and_formats import (
    SAMDirFmt, SAMFormat, _SKIP_FLAGS)

# samtools filter expressions. SAM uses '*' where SEQ or QUAL is unavailable,
# which samtools reads as an empty sequence, and as quality scores of 255
# that no real quality score reaches.
_HAS_SEQUENCE = 'length(seq) > 0'
_MISSING_QUALITY = '%s && min(qual) == 255' % _HAS_SEQUENCE


# A SAMDirFmt holds the alignments of a single sample, one file per
# reference, so each reference's reads become one entry of the per-sample
# sequence formats, identified by the reference's ID.
@plugin.register_transformer
def _1(dirfmt: SAMDirFmt) -> SingleLanePerSampleSingleEndFastqDirFmt:
    result = SingleLanePerSampleSingleEndFastqDirFmt()
    manifest_data = []

    for index, (reference, sam) in enumerate(_iter_references(dirfmt)):
        # samtools would otherwise fill in a default quality score for these
        # reads, rather than refusing to write them.
        missing = int(_samtools_view(sam, _MISSING_QUALITY, '-c').stdout)
        if missing:
            raise ValueError(
                '%d read(s) aligned to %r have no quality scores, so they '
                'cannot be represented as SampleData[SequencesWithQuality]. '
                'View this data as SampleData[Sequences] instead.'
                % (missing, reference))

        output_fp = result.sequences.path_maker(
            sample_id=reference,
            # These are not used internally by QIIME 2, so their values
            # don't matter beyond keeping the filenames unique.
            barcode_id=index, lane_number=1, read_number=1)

        # samtools compresses the output because its name ends in '.gz'.
        _samtools_view(sam, _HAS_SEQUENCE,
                       '--output-fmt', 'fastq', '-o', str(output_fp))

        manifest_data.append([reference, output_fp.name, 'forward'])

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
        for reference, sam in _iter_references(dirfmt):
            fasta = _samtools_view(
                sam, _HAS_SEQUENCE, '--output-fmt', 'fasta').stdout

            # samtools writes every sequence on the line after its header.
            # QIIME 1 demultiplexed FASTA identifies a sequence by
            # '<sample-id>_<seq-id>', where seq-id is unique within the
            # sample, so the read name samtools writes is replaced.
            for seq_id, sequence in enumerate(fasta.splitlines()[1::2]):
                fh.write('>%s_%d\n%s\n' % (reference, seq_id, sequence))

    return result


def _iter_references(dirfmt):
    for path, sam in sorted(dirfmt.sams.iter_views(SAMFormat)):
        yield path.stem, str(sam)


def _samtools_view(sam_path, expression, *options):
    # Only primary alignments are kept, so that every read is seen once. When
    # writing FASTQ or FASTA, samtools also undoes the reverse complement
    # that SAM stores for reads that aligned to the reverse strand.
    return _run_command([
        'samtools', 'view', '-F', str(_SKIP_FLAGS), '-e', expression,
        *options, sam_path])
