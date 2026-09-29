# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import qiime2
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import PairedEndSequencesWithQuality

from q2_multi_amplicon._methods import _classify_reads_defaults

# These mirror the defaults of the q2-cutadapt and q2-vsearch actions that
# this pipeline calls, so that `qiime multi-amplicon multi-amplicon-analysis`
# reports the values that will actually be used.
_trim_defaults = {
    'metadata': None,
    'adapter_f': None,
    'adapter_r': None,
    'front_f': None,
    'front_r': None,
    'anywhere_f': None,
    'anywhere_r': None,
    'forward_cut': 0,
    'reverse_cut': 0,
    'error_rate': 0.1,
    'indels': True,
    'times': 1,
    'overlap': 3,
    'match_read_wildcards': False,
    'match_adapter_wildcards': True,
    'minimum_length': 1,
    'discard_untrimmed': False,
    'quality_cutoff_5end': 0,
    'quality_cutoff_3end': 0,
}

# trim_single names its adapter parameters, and the metadata columns that
# can stand in for them, without a read direction. This pipeline uses the
# paired-end names for both kinds of input, so the forward-read columns are
# renamed before single-end reads are trimmed. Like cutadapt, a column name
# may use '-' in place of '_'.
_single_end_adapter_names = {
    'adapter_f': 'adapter',
    'front_f': 'front',
    'anywhere_f': 'anywhere',
}
_reverse_read_adapter_names = ('adapter_r', 'front_r', 'anywhere_r')

_merge_defaults = {
    'truncqual': None,
    'minlen': 1,
    'maxns': None,
    'allowmergestagger': False,
    'minovlen': 10,
    'maxdiffs': 10,
    'minmergelen': None,
    'maxmergelen': None,
    'maxee': None,
}

_call_alleles_defaults = {
    'reference_alleles': None,
    'min_coverage': 10,
    'proportion': 0.5,
}


def multi_amplicon_analysis(
        ctx, sequences, database,
        reference_alleles=_call_alleles_defaults['reference_alleles'],
        metadata=_trim_defaults['metadata'],
        adapter_f=_trim_defaults['adapter_f'],
        front_f=_trim_defaults['front_f'],
        anywhere_f=_trim_defaults['anywhere_f'],
        adapter_r=_trim_defaults['adapter_r'],
        front_r=_trim_defaults['front_r'],
        anywhere_r=_trim_defaults['anywhere_r'],
        forward_cut=_trim_defaults['forward_cut'],
        reverse_cut=_trim_defaults['reverse_cut'],
        error_rate=_trim_defaults['error_rate'],
        indels=_trim_defaults['indels'],
        times=_trim_defaults['times'],
        overlap=_trim_defaults['overlap'],
        match_read_wildcards=_trim_defaults['match_read_wildcards'],
        match_adapter_wildcards=_trim_defaults['match_adapter_wildcards'],
        minimum_length=_trim_defaults['minimum_length'],
        discard_untrimmed=_trim_defaults['discard_untrimmed'],
        quality_cutoff_5end=_trim_defaults['quality_cutoff_5end'],
        quality_cutoff_3end=_trim_defaults['quality_cutoff_3end'],
        truncqual=_merge_defaults['truncqual'],
        minlen=_merge_defaults['minlen'],
        maxns=_merge_defaults['maxns'],
        allowmergestagger=_merge_defaults['allowmergestagger'],
        minovlen=_merge_defaults['minovlen'],
        maxdiffs=_merge_defaults['maxdiffs'],
        minmergelen=_merge_defaults['minmergelen'],
        maxmergelen=_merge_defaults['maxmergelen'],
        maxee=_merge_defaults['maxee'],
        mode=_classify_reads_defaults['mode'],
        sensitivity=_classify_reads_defaults['sensitivity'],
        min_coverage=_call_alleles_defaults['min_coverage'],
        proportion=_call_alleles_defaults['proportion'],
        n_threads=_classify_reads_defaults['n_threads']):
    paired = sequences.type <= SampleData[PairedEndSequencesWithQuality]

    shared_trim_params = {
        'error_rate': error_rate,
        'indels': indels,
        'times': times,
        'overlap': overlap,
        'match_read_wildcards': match_read_wildcards,
        'match_adapter_wildcards': match_adapter_wildcards,
        'minimum_length': minimum_length,
        'discard_untrimmed': discard_untrimmed,
        'quality_cutoff_5end': quality_cutoff_5end,
        'quality_cutoff_3end': quality_cutoff_3end,
        'cores': n_threads,
    }

    if paired:
        trim_paired = ctx.get_action('cutadapt', 'trim_paired')
        merge_pairs = ctx.get_action('vsearch', 'merge_pairs')

        trimmed, trim_stats = trim_paired(
            demultiplexed_sequences=sequences, metadata=metadata,
            adapter_f=adapter_f, front_f=front_f, anywhere_f=anywhere_f,
            adapter_r=adapter_r, front_r=front_r, anywhere_r=anywhere_r,
            forward_cut=forward_cut, reverse_cut=reverse_cut,
            **shared_trim_params)

        # merge_pairs also reports the pairs it could not merge. Those reads
        # cannot be classified against the index, so they are dropped here
        # rather than returned as an output that would be empty for
        # single-end input.
        classified_sequences, _ = merge_pairs(
            demultiplexed_seqs=trimmed, truncqual=truncqual, minlen=minlen,
            maxns=maxns, allowmergestagger=allowmergestagger,
            minovlen=minovlen, maxdiffs=maxdiffs, minmergelen=minmergelen,
            maxmergelen=maxmergelen, maxee=maxee, threads=n_threads)
    else:
        _reject_paired_only_params(
            metadata=metadata, adapter_r=adapter_r, front_r=front_r,
            anywhere_r=anywhere_r, reverse_cut=reverse_cut,
            truncqual=truncqual, minlen=minlen,
            maxns=maxns, allowmergestagger=allowmergestagger,
            minovlen=minovlen, maxdiffs=maxdiffs, minmergelen=minmergelen,
            maxmergelen=maxmergelen, maxee=maxee)

        trim_single = ctx.get_action('cutadapt', 'trim_single')
        classified_sequences, trim_stats = trim_single(
            demultiplexed_sequences=sequences,
            metadata=_as_single_end_metadata(metadata),
            adapter=adapter_f, front=front_f, anywhere=anywhere_f,
            cut=forward_cut, **shared_trim_params)

    classify_reads = ctx.get_action('multi_amplicon', 'classify_reads')
    classifications, = classify_reads(
        sequences=classified_sequences, database=database,
        n_threads=n_threads, mode=mode, sensitivity=sensitivity)

    count_classifications = ctx.get_action(
        'multi_amplicon', 'count_classifications')
    table, = count_classifications(classifications=classifications)

    dereplicate_classifications = ctx.get_action(
        'multi_amplicon', 'dereplicate_classifications')
    dereplicated_tables, dereplicated_sequences = \
        dereplicate_classifications(classifications=classifications)

    call_alleles = ctx.get_action('real_allele', 'call_alleles')
    genotypes = {}
    for locus, locus_table in dereplicated_tables.items():
        # call_alleles fails when no sequence has min_coverage reads. With
        # many amplicons sequenced together some are expected to be that
        # sparse, so they are left out instead of failing the whole run.
        counts = locus_table.view(pd.DataFrame)
        if not (counts >= min_coverage).values.any():
            continue

        # Alleles are named per locus, so each locus is matched with its own
        # reference alleles, if any were given for it.
        locus_reference = None
        if (reference_alleles is not None
                and locus in reference_alleles.keys()):
            locus_reference = reference_alleles[locus]

        genotypes[locus] = call_alleles(
            sequences=dereplicated_sequences[locus], table=locus_table,
            reference_alleles=locus_reference, min_coverage=min_coverage,
            proportion=proportion)

    def collect(output):
        return {locus: getattr(results, output)
                for locus, results in genotypes.items()}

    return (classified_sequences, classifications, table,
            dereplicated_tables, dereplicated_sequences,
            collect('alleles'), collect('genotyped_sequences'),
            collect('genotype_table'), collect('allele_frequencies'),
            collect('heterozygosity'), trim_stats)


def _reject_paired_only_params(metadata, **provided):
    # The reverse-read and merging parameters mean nothing when the reads are
    # single-end. If any of them was set the user has misunderstood what the
    # pipeline is about to do, so fail rather than silently ignore them.
    defaults = dict(_trim_defaults)
    defaults.update(_merge_defaults)

    set_params = [p for p, v in provided.items() if v != defaults[p]]
    if metadata is not None:
        set_params += [
            '%s (metadata column)' % c for c in metadata.columns
            if c.replace('-', '_') in _reverse_read_adapter_names]
    set_params.sort()

    if set_params:
        raise ValueError(
            'The following parameter(s) only apply to paired-end sequences, '
            'but single-end sequences were provided: %s. Remove these '
            'parameters, or provide SampleData[PairedEndSequencesWithQuality] '
            'as input.' % ', '.join(set_params))


def _as_single_end_metadata(metadata):
    if metadata is None:
        return None

    renames = {
        c: _single_end_adapter_names[c.replace('-', '_')]
        for c in metadata.columns
        if c.replace('-', '_') in _single_end_adapter_names}
    if not renames:
        return metadata

    return qiime2.Metadata(metadata.to_dataframe().rename(columns=renames))
