# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.plugin import (Bool, Choices, Citations, Collection, Float, Int,
                           List, Metadata, Plugin, Range, Str, Threads,
                           TypeMap)
from q2_types.bowtie2 import Bowtie2Index
from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.metadata import ImmutableMetadata
from q2_types.per_sample_sequences import (
    JoinedSequencesWithQuality, PairedEndSequencesWithQuality,
    Sequences, SequencesWithQuality)
from q2_types.sample_data import SampleData

from q2_multi_amplicon import __version__
from q2_multi_amplicon._methods import (
    classify_reads, count_classifications, dereplicate_classifications,
    extend_index)
from q2_multi_amplicon._pipelines import multi_amplicon_analysis
from q2_multi_amplicon._types_and_formats import (
    SAMDirFmt, SAMFormat, SAMOutput)

citations = Citations.load("citations.bib", package="q2_multi_amplicon")

plugin = Plugin(
    name="multi-amplicon",
    version=__version__,
    website="https://example.com/q2-multi-amplicon",
    package="q2_multi_amplicon",
    description=("A QIIME 2 plugin for working with multiple amplicon "
                 "regions/datasets."),
    short_description="Multi-amplicon workflows.",
    # The plugin-level citation of 'Caporaso-Bolyen-2024' is provided as
    # an example. You can replace this with citations to other references
    # in citations.bib.
    citations=[citations['Caporaso-Bolyen-2024']]
)

plugin.register_formats(SAMFormat, SAMDirFmt)
plugin.register_semantic_types(SAMOutput)
plugin.register_semantic_type_to_format(
    SampleData[SAMOutput], artifact_format=SAMDirFmt)

_classify_reads_parameters = {
    'n_threads': Threads,
    'mode': Str % Choices(['local', 'global']),
    'sensitivity': Str % Choices(
        ['very-fast', 'fast', 'sensitive', 'very-sensitive']),
}

_classify_reads_parameter_descriptions = {
    'n_threads': 'Number of alignment threads to launch.',
    'mode': 'Bowtie2 alignment settings. See the bowtie2 manual for more '
            'details.',
    'sensitivity': 'Bowtie2 alignment sensitivity. See the bowtie2 manual '
                   'for more details.',
}

plugin.methods.register_function(
    function=classify_reads,
    inputs={
        'sequences': SampleData[SequencesWithQuality |
                                JoinedSequencesWithQuality],
        'database': Bowtie2Index,
    },
    parameters=_classify_reads_parameters,
    outputs=[('classifications', SampleData[SAMOutput])],
    input_descriptions={
        'sequences': 'The single-end (or merged) sequences to classify.',
        'database': 'Bowtie2 index of the reference sequences that reads '
                    'will be assigned to.',
    },
    parameter_descriptions=_classify_reads_parameter_descriptions,
    output_descriptions={
        'classifications': 'One SAM file per sample, recording the reference '
                           'sequence each read was assigned to.',
    },
    name='Classify reads against a Bowtie2 index',
    description=(
        'Align each sample\'s reads to a Bowtie2 index with bowtie2, and '
        'return the resulting alignments as one SAM file per sample. Reads '
        'that did not align to any reference are left out of the output. '
        'The alignments can be viewed as SampleData[Sequences] or '
        'SampleData[SequencesWithQuality] to recover the classified reads '
        'themselves.'),
    citations=[citations['langmead2012fast']]
)

plugin.methods.register_function(
    function=count_classifications,
    inputs={'classifications': SampleData[SAMOutput]},
    parameters={},
    outputs=[('table', FeatureTable[Frequency])],
    input_descriptions={
        'classifications': 'The per-sample alignments to count.',
    },
    parameter_descriptions={},
    output_descriptions={
        'table': 'The number of reads in each sample that were assigned to '
                 'each reference sequence.',
    },
    name='Count the reads assigned to each reference sequence',
    description=(
        'Count, per sample, the reads that were assigned to each reference '
        'sequence, producing a feature table whose features are the IDs of '
        'the reference sequences. Only each read\'s primary alignment is '
        'counted, so a read is counted at most once. Every sample is '
        'included in the table, even one in which no read was assigned to '
        'any reference.'),
    citations=[]
)

_dereplicated_output_descriptions = {
    'dereplicated_tables': 'One feature table per reference sequence, keyed '
                           'by the reference\'s ID, counting the reads of '
                           'each sample that were assigned to that '
                           'reference, per distinct read sequence.',
    'dereplicated_sequences': 'The distinct read sequences assigned to each '
                              'reference sequence, keyed by the reference\'s '
                              'ID and oriented like the reference. Each '
                              'feature ID is the MD5 hash of its sequence, '
                              'matching the table with the same key.',
}

plugin.methods.register_function(
    function=dereplicate_classifications,
    inputs={'classifications': SampleData[SAMOutput]},
    parameters={},
    outputs=[
        ('dereplicated_tables', Collection[FeatureTable[Frequency]]),
        ('dereplicated_sequences', Collection[FeatureData[Sequence]]),
    ],
    input_descriptions={
        'classifications': 'The per-sample alignments to dereplicate.',
    },
    parameter_descriptions={},
    output_descriptions=_dereplicated_output_descriptions,
    name='Dereplicate the reads assigned to each reference sequence',
    description=(
        'Group each sample\'s reads by the reference sequence they were '
        'assigned to, and collapse identical reads within each group into '
        'one feature. Each reference gets its own feature table and '
        'sequences, which together give the number of reads with each '
        'distinct sequence in each sample, as used for allele calling. '
        'Reads that aligned to the reverse strand are reverse complemented '
        'in SAM, so all of a reference\'s sequences share its orientation. '
        'Only each read\'s primary alignment is used. References that no '
        'read was assigned to are left out.'),
    citations=[]
)

plugin.methods.register_function(
    function=extend_index,
    inputs={
        'sequences': SampleData[SequencesWithQuality | Sequences],
        'database': Bowtie2Index,
    },
    parameters={'n_threads': Threads},
    outputs=[('extended_database', Bowtie2Index)],
    input_descriptions={
        'sequences': 'The sequences to add to the index.',
        'database': 'The Bowtie2 index to extend.',
    },
    parameter_descriptions={
        'n_threads': 'Number of threads to build the new index with.',
    },
    output_descriptions={
        'extended_database': 'A new Bowtie2 index built from the sequences '
                             'of the original index plus the provided '
                             'sequences.',
    },
    name='Add sequences to a Bowtie2 index',
    description=(
        'Write the reference sequences of an existing Bowtie2 index back out '
        'as FASTA with bowtie2-inspect, append the provided sequences to '
        'them, and build a new index from the combination. The original '
        'index is left untouched. Adding a sequence whose ID is already in '
        'the index is an error, because the resulting index could not '
        'distinguish the two.'),
    citations=[citations['langmead2012fast']]
)

_multi_amplicon_analysis_parameters = {
    'metadata': Metadata,
    'adapter_f': List[Str],
    'front_f': List[Str],
    'anywhere_f': List[Str],
    'adapter_r': List[Str],
    'front_r': List[Str],
    'anywhere_r': List[Str],
    'forward_cut': Int,
    'reverse_cut': Int,
    'error_rate': Float % Range(0, 1, inclusive_start=True,
                                inclusive_end=True),
    'indels': Bool,
    'times': Int % Range(1, None),
    'overlap': Int % Range(1, None),
    'match_read_wildcards': Bool,
    'match_adapter_wildcards': Bool,
    'minimum_length': Int % Range(1, None),
    'discard_untrimmed': Bool,
    'quality_cutoff_5end': Int % Range(0, None),
    'quality_cutoff_3end': Int % Range(0, None),
    'truncqual': Int % Range(0, None),
    'minlen': Int % Range(0, None),
    'maxns': Int % Range(0, None),
    'allowmergestagger': Bool,
    'minovlen': Int % Range(5, None),
    'maxdiffs': Int % Range(0, None),
    'minmergelen': Int % Range(0, None),
    'maxmergelen': Int % Range(0, None),
    'maxee': Float % Range(0., None),
    'min_coverage': Int % Range(1, None),
    'proportion': Float % Range(0, 1, inclusive_end=True),
}
_multi_amplicon_analysis_parameters.update(_classify_reads_parameters)

_multi_amplicon_analysis_parameter_descriptions = {
    'metadata': 'Adapters to trim, given as metadata instead of as '
                'parameters. Each column is named after one of the adapter '
                'parameters (`adapter_f`, `front_f`, `anywhere_f`, '
                '`adapter_r`, `front_r`, `anywhere_r`) and holds one adapter '
                'per row. Columns may be combined with the adapter '
                'parameters, but the same adapter parameter may not be set '
                'both ways. Reverse-read columns are only allowed for '
                'paired-end input.',
    'adapter_f': 'Sequence of an adapter ligated to the 3\' end of the '
                 'forward read. Passed to cutadapt as `adapter` for '
                 'single-end input.',
    'front_f': 'Sequence of an adapter ligated to the 5\' end of the forward '
               'read. Passed to cutadapt as `front` for single-end input.',
    'anywhere_f': 'Sequence of an adapter that may be ligated to either end '
                  'of the forward read. Passed to cutadapt as `anywhere` for '
                  'single-end input.',
    'adapter_r': 'Sequence of an adapter ligated to the 3\' end of the '
                 'reverse read. Paired-end input only.',
    'front_r': 'Sequence of an adapter ligated to the 5\' end of the reverse '
               'read. Paired-end input only.',
    'anywhere_r': 'Sequence of an adapter that may be ligated to either end '
                  'of the reverse read. Paired-end input only.',
    'forward_cut': 'Unconditionally remove bases from the forward read '
                   'before adapter trimming. Positive values trim from the '
                   '5\' end, negative values from the 3\' end.',
    'reverse_cut': 'Unconditionally remove bases from the reverse read '
                   'before adapter trimming. Paired-end input only.',
    'error_rate': 'Maximum allowed error rate when matching adapters.',
    'indels': 'Allow insertions or deletions of bases when matching '
              'adapters.',
    'times': 'Remove multiple occurrences of an adapter if it is repeated, '
             'up to this many times.',
    'overlap': 'Require at least this many bases of overlap between the read '
               'and the adapter for an adapter to be found.',
    'match_read_wildcards': 'Interpret IUPAC wildcards in the reads.',
    'match_adapter_wildcards': 'Interpret IUPAC wildcards in the adapters.',
    'minimum_length': 'Discard reads shorter than this after trimming. '
                      'Reads of length zero are always discarded, because '
                      'downstream actions cannot handle them.',
    'discard_untrimmed': 'Discard reads in which no adapter was found.',
    'quality_cutoff_5end': 'Trim low-quality bases from the 5\' end of each '
                           'read before adapter removal.',
    'quality_cutoff_3end': 'Trim low-quality bases from the 3\' end of each '
                           'read before adapter removal.',
    'truncqual': 'Truncate reads at the first base with this quality score '
                 'or lower before merging. Paired-end input only.',
    'minlen': 'Reads shorter than this after truncation are discarded before '
              'merging. Paired-end input only.',
    'maxns': 'Reads with more than this many N characters are discarded '
             'before merging. Paired-end input only.',
    'allowmergestagger': 'Allow merging of staggered read pairs. Paired-end '
                         'input only.',
    'minovlen': 'Minimum length of the area of overlap between reads during '
                'merging. Paired-end input only.',
    'maxdiffs': 'Maximum number of mismatches in the area of overlap during '
                'merging. Paired-end input only.',
    'minmergelen': 'Minimum length of a merged read for it to be retained. '
                   'Paired-end input only.',
    'maxmergelen': 'Maximum length of a merged read for it to be retained. '
                   'Paired-end input only.',
    'maxee': 'Maximum number of expected errors in a merged read for it to '
             'be retained. Paired-end input only.',
    'n_threads': 'Number of threads to use. Passed on to cutadapt, vsearch '
                 'and bowtie2.',
    'mode': _classify_reads_parameter_descriptions['mode'],
    'sensitivity': _classify_reads_parameter_descriptions['sensitivity'],
    'min_coverage': 'Distinct read sequences supported by fewer reads than '
                    'this in a sample are not considered when calling that '
                    'sample\'s alleles. Reference sequences with no '
                    'sequence this well supported in any sample are left '
                    'out of the allele outputs.',
    'proportion': 'Only the two most abundant sequences of each sample at a '
                  'reference sequence are considered. The second is kept, '
                  'making the sample heterozygous there, if its read count '
                  'is at least this proportion of the first\'s; otherwise the '
                  'sample is homozygous for the most abundant sequence.',
}

# Single-end reads come back out single-end, while paired-end reads are
# merged along the way and so come back out joined.
_input_seqs, _output_seqs = TypeMap({
    SampleData[SequencesWithQuality]: SampleData[SequencesWithQuality],
    SampleData[PairedEndSequencesWithQuality]:
        SampleData[JoinedSequencesWithQuality],
})

plugin.pipelines.register_function(
    function=multi_amplicon_analysis,
    inputs={
        'sequences': _input_seqs,
        'database': Bowtie2Index,
        'reference_alleles': Collection[FeatureData[Sequence]],
    },
    parameters=_multi_amplicon_analysis_parameters,
    outputs=[
        ('classified_sequences', _output_seqs),
        ('classifications', SampleData[SAMOutput]),
        ('table', FeatureTable[Frequency]),
        ('dereplicated_tables', Collection[FeatureTable[Frequency]]),
        ('dereplicated_sequences', Collection[FeatureData[Sequence]]),
        ('alleles', Collection[FeatureData[Sequence]]),
        ('genotyped_sequences', Collection[FeatureData[Sequence]]),
        ('genotype_tables', Collection[ImmutableMetadata]),
        ('allele_frequencies', Collection[ImmutableMetadata]),
        ('heterozygosity', Collection[ImmutableMetadata]),
        ('trim_stats', ImmutableMetadata),
    ],
    input_descriptions={
        'sequences': 'The demultiplexed single- or paired-end sequences to '
                     'analyze.',
        'database': 'Bowtie2 index of the reference sequences that reads '
                    'will be assigned to. Build one with '
                    '`qiime quality-control bowtie2-build`.',
        'reference_alleles': 'Previously named alleles of each reference '
                             'sequence, keyed by the reference\'s ID, such '
                             'as the `alleles` output of an earlier run. A '
                             'sequence identical to a reference allele is '
                             'given its name, and new alleles are numbered '
                             'after the highest reference allele. Reference '
                             'sequences without an entry get new allele '
                             'names starting from A1.',
    },
    parameter_descriptions=_multi_amplicon_analysis_parameter_descriptions,
    output_descriptions={
        'classified_sequences': 'The trimmed reads (merged, if the input was '
                                'paired-end) that were classified.',
        'classifications': 'One SAM file per sample, recording the reference '
                           'sequence each read was assigned to.',
        'table': 'The number of reads in each sample that were assigned to '
                 'each reference sequence.',
        **_dereplicated_output_descriptions,
        'alleles': 'The allele sequences of each reference sequence, keyed by '
                   'the reference\'s ID and named A1, A2, ... in order of '
                   'decreasing frequency. Where reference alleles were '
                   'given, these contain all of them followed by any new '
                   'alleles, so they can be the reference alleles of the '
                   'next run.',
        'genotyped_sequences': 'The sequences that the genotypes of each '
                               'reference sequence were called from, keyed '
                               'by the reference\'s ID, with IDs of the form '
                               '<sample>_<read count>_<allele>.',
        'genotype_tables': 'The genotype of each sample at each reference '
                           'sequence, keyed by the reference\'s ID: the read '
                           'count and allele of its most abundant (1) and '
                           'second (2) sequence. Homozygous samples list the '
                           'same read count and allele twice.',
        'allele_frequencies': 'The number of genotyped sequences carrying '
                              'each allele, per reference sequence.',
        'heterozygosity': 'The number and percentage of samples that are '
                          'homozygous and heterozygous, per reference '
                          'sequence.',
        'trim_stats': 'Per-sample summary of the cutadapt trimming step.',
    },
    name='Multi-amplicon analysis',
    description=(
        'Trim adapters and primers from demultiplexed single- or paired-end '
        'reads with cutadapt, merge read pairs with vsearch if the input was '
        'paired-end, classify the resulting reads per sample against a '
        'Bowtie2 index of reference sequences, and count the reads assigned '
        'to each reference sequence in each sample, both in total and per '
        'distinct read sequence. Finally, call alleles and genotypes at each '
        'reference sequence from those distinct sequences with '
        '`qiime real-allele call-alleles`. Pairs that vsearch '
        'could not merge are not classified and are not returned. '
        'Reverse-read and merging parameters, and reverse-read metadata '
        'columns, may only be provided when the input is paired-end.'),
    citations=[]
)

importlib.import_module('q2_multi_amplicon._transformers')
