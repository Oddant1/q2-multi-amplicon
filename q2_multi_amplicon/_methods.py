# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import hashlib
import itertools
import shutil
import subprocess
import tempfile
from pathlib import Path

import biom
import numpy as np
import pandas as pd
from qiime2.plugin import get_available_cores
from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt, QIIME1DemuxDirFmt,
    QIIME1DemuxFormat)

from q2_multi_amplicon._types_and_formats import (
    SAMDirFmt, SAMFormat, _UNAVAILABLE, _iter_primary_records)


_classify_reads_defaults = {
    'n_threads': 1,
    'mode': 'local',
    'sensitivity': 'sensitive',
}

_extend_index_defaults = {
    'n_threads': 1,
}

# bowtie2-build names every file in the index after this stem. It matches
# what `qiime quality-control bowtie2-build` uses.
_INDEX_BASENAME = 'db'


def classify_reads(
        sequences: CasavaOneEightSingleLanePerSampleDirFmt,
        database: Bowtie2IndexDirFmt,
        n_threads: int = _classify_reads_defaults['n_threads'],
        mode: str = _classify_reads_defaults['mode'],
        sensitivity: str = _classify_reads_defaults['sensitivity'],
) -> SAMDirFmt:
    if n_threads == 0:
        n_threads = get_available_cores()

    index_base = str(database.path / database.get_basename())
    preset = _bowtie2_preset(mode, sensitivity)

    classifications = {}
    for sample_id, forward, _ in sequences.manifest.itertuples():
        with tempfile.TemporaryDirectory(
                prefix='q2-multi-amplicon-') as tmpdir:
            unsorted_sam = str(Path(tmpdir) / 'unsorted.sam')
            sorted_sam = str(Path(tmpdir) / 'sorted.sam')

            # --no-unal keeps reads that could not be assigned to any
            # reference out of the output entirely.
            _run_command([
                'bowtie2', '-p', str(n_threads), preset, '--no-unal',
                '-x', index_base, '-U', str(forward), '-S', unsorted_sam])

            # Sorting brings each reference's records together, so that they
            # can be split into one file per reference while holding only
            # one of those files open at a time.
            _run_command([
                'samtools', 'sort', '-O', 'sam', '-T',
                str(Path(tmpdir) / 'sort'), '-o', sorted_sam, unsorted_sam])

            # A sample none of whose reads were assigned is still included,
            # with no files, so that it is accounted for downstream.
            classifications[sample_id] = _split_by_reference(sorted_sam)

    return classifications


def count_classifications(classifications: SAMDirFmt) -> biom.Table:
    counts = {}
    for sample_id, alignments in sorted(classifications.items()):
        # A sample none of whose reads were assigned still gets a column, so
        # that every sample that was classified is accounted for.
        counts[sample_id] = collections.Counter(
            fields[2] for fields in _iter_sample_records(alignments))

    table = pd.DataFrame(counts).fillna(0).astype(int).sort_index()
    return biom.Table(table.values, list(table.index), list(table.columns))


def dereplicate_classifications(
        classifications: SAMDirFmt) -> (biom.Table, pd.Series):
    sample_ids = []
    # reference -> read sequence -> sample -> number of reads
    counts = collections.defaultdict(
        lambda: collections.defaultdict(collections.Counter))
    for sample_id, alignments in sorted(classifications.items()):
        sample_ids.append(sample_id)
        for fields in _iter_sample_records(alignments):
            # SAM stores a read that aligned to the reverse strand as its
            # reverse complement, so every sequence here is already in the
            # orientation of its reference, whichever strand it came from.
            reference, sequence = fields[2], fields[9]
            if sequence != _UNAVAILABLE:
                counts[reference][sequence][sample_id] += 1

    tables, sequences = {}, {}
    for reference, by_sequence in sorted(counts.items()):
        # Most abundant first, with ties broken by sequence so that the
        # order does not depend on the order of the reads.
        ordered = sorted(by_sequence,
                         key=lambda s: (-sum(by_sequence[s].values()), s))
        ids = [hashlib.md5(s.encode()).hexdigest() for s in ordered]

        data = np.array([[by_sequence[s][sample_id]
                          for sample_id in sample_ids] for s in ordered])
        tables[reference] = biom.Table(data, ids, list(sample_ids))
        sequences[reference] = pd.Series(ordered, index=ids)

    return tables, sequences


def extend_index(
        sequences: QIIME1DemuxDirFmt,
        database: Bowtie2IndexDirFmt,
        n_threads: int = _extend_index_defaults['n_threads'],
) -> Bowtie2IndexDirFmt:
    if n_threads == 0:
        n_threads = get_available_cores()

    index_base = str(database.path / database.get_basename())
    new_sequences = str(sequences.file.view(QIIME1DemuxFormat))

    with tempfile.TemporaryDirectory(prefix='q2-multi-amplicon-') as tmpdir:
        combined = Path(tmpdir) / 'combined.fasta'

        # bowtie2-inspect writing the indexed references back out as FASTA is
        # the only way to recover the sequences an index was built from.
        with open(combined, 'w') as fh:
            _run_command(['bowtie2-inspect', index_base], stdout=fh)

        _reject_duplicate_ids(combined, new_sequences)

        with open(combined, 'a') as out_fh, open(new_sequences) as in_fh:
            shutil.copyfileobj(in_fh, out_fh)

        new_database = Bowtie2IndexDirFmt()
        _run_command([
            'bowtie2-build', '--threads', str(n_threads), str(combined),
            str(new_database.path / _INDEX_BASENAME)])

    return new_database


def _reject_duplicate_ids(exported_path, new_sequences_path):
    # Two references sharing a name in one index would make any assignment to
    # that name ambiguous, and bowtie2-build will not complain about it.
    duplicates = sorted(
        _fasta_ids(exported_path) & _fasta_ids(new_sequences_path))

    if duplicates:
        raise ValueError(
            'The following sequence ID(s) are already in the index, so '
            'adding them would make the resulting index ambiguous: %s. '
            'Rename these sequences, or remove them from the input.'
            % ', '.join(duplicates))


def _fasta_ids(path):
    ids = set()
    with open(str(path)) as fh:
        for line in fh:
            if line.startswith('>'):
                fields = line[1:].split()
                if fields:
                    ids.add(fields[0])
    return ids


def _bowtie2_preset(mode, sensitivity):
    # bowtie2 spells its local presets '--sensitive-local' and its end-to-end
    # presets '--sensitive', so the mode is only part of the flag in local
    # mode.
    if mode == 'local':
        return '--%s-%s' % (sensitivity, mode)
    return '--%s' % sensitivity


def _split_by_reference(sam_path):
    before, sq_lines, after = _read_header(sam_path)

    alignments = SAMDirFmt()
    with open(sam_path) as fh:
        records = (line for line in fh if not line.startswith('@'))
        # The records are sorted, so each reference's records are
        # contiguous and each reference's file is written in one go.
        for reference, lines in itertools.groupby(records, key=_rname):
            path = alignments.sams.path_maker(reference=reference)
            with open(path, 'w') as out_fh:
                # Only the reference's own @SQ line is kept, so that the file
                # describes only the reference whose alignments it holds.
                out_fh.writelines(before + [sq_lines[reference]] + after)
                out_fh.writelines(lines)

    return alignments


def _read_header(sam_path):
    """Split a SAM header into its @SQ lines, keyed by reference, and the
    header lines before and after them."""
    before, sq_lines, after = [], {}, []
    with open(sam_path) as fh:
        for line in fh:
            if not line.startswith('@'):
                break
            if line.startswith('@SQ\t'):
                sq_lines[_sq_name(line)] = line
            elif sq_lines:
                after.append(line)
            else:
                before.append(line)
    return before, sq_lines, after


def _sq_name(line):
    for field in line.rstrip('\n').split('\t')[1:]:
        if field.startswith('SN:'):
            return field[3:]
    raise ValueError('SAM @SQ header line has no SN field: %r' % line)


def _rname(line):
    return line.split('\t', 3)[2]


def _iter_sample_records(alignments):
    for _, sam in sorted(alignments.sams.iter_views(SAMFormat)):
        yield from _iter_primary_records(str(sam))


def _run_command(cmd, stdout=None):
    try:
        return subprocess.run(cmd, check=True, text=True,
                              stdout=stdout if stdout is not None
                              else subprocess.PIPE,
                              stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        detail = (e.stderr or '').strip() or (e.stdout or '').strip()
        raise RuntimeError(
            '%s failed with exit code %d: %s'
            % (cmd[0], e.returncode, detail)) from e
