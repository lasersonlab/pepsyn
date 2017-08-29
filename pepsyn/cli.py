# Copyright 2016 Uri Laserson
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from logging import captureWarnings

# biopython has a bunch of annoying warnings bc Seq comparisons changed
captureWarnings(True)

from click import group, command, option, argument, File, Choice
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Data.CodonTable import standard_dna_table

from pepsyn.operations import (
    reverse_translate, recode_site_from_cds, x_to_ggsg, disambiguate_iupac_aa,
    tile as tile_op, ctermpep as cterm_oligo, pad_ggsg)
from pepsyn.codons import (
    FreqWeightedCodonSampler, UniformCodonSampler, ecoli_codon_usage,
    zero_non_amber_stops, zero_low_freq_codons)
from pepsyn.util import site2dna


@group(context_settings={'help_option_names': ['-h', '--help']})
def cli():
    """pepsyn -- peptide synthesis design"""
    pass


# reusable args
argument_input = argument('input', type=File('r'))
argument_output = argument('output', type=File('w'))


@cli.command()
@argument_input
@argument_output
@option('--length', '-l', type=int, help='Length of output oligos')
@option('--overlap', '-p', type=int, help='Overlap of oligos')
def tile(input, output, length, overlap):
    """tile a set of sequences"""
    for seqrecord in SeqIO.parse(input, 'fasta'):
        for (start, end, t) in tile_op(seqrecord.seq, length, overlap):
                output_title = '{}|{}-{}'.format(seqrecord.id, start, end)
                output_record = SeqRecord(t, output_title, description='')
                SeqIO.write(output_record, output, 'fasta')


@cli.command()
@argument_input
@argument_output
@option('--length', '-l', type=int, help='Length of output C-terminal oligo')
@option('--add-stop', '-s', is_flag=True, help='Add a stop codon to peptide')
def ctermpep(input, output, length, add_stop):
    """extract C-terminal peptide (AA alphabets)"""
    for seqrecord in SeqIO.parse(input, 'fasta'):
        oligo = cterm_oligo(seqrecord.seq, length, add_stop=add_stop)
        output_title = '{}|CTERM'.format(seqrecord.id)
        if add_stop:
            output_title = '{}|STOP'.format(output_title)
        output_record = SeqRecord(oligo, output_title, description='')
        SeqIO.write(output_record, output, 'fasta')


@cli.command()
@argument_input
@argument_output
def stripstop(input, output):
    """strip stop codons from end of protein sequence"""
    for seqrecord in SeqIO.parse(input, 'fasta'):
        seqrecord.seq = seqrecord.seq.rstrip('*')
        SeqIO.write(seqrecord, output, 'fasta')


@cli.command()
@argument_input
@argument_output
@option('--length', '-l', type=int, help='Target length for peptide')
@option('--n-term', '-n', 'terminus', flag_value='N',
        help='Pad the N-terminus')
@option('--c-term', '-c', 'terminus', flag_value='C', default=True,
        help='Pad the C-terminus')
def pad(input, output, length, terminus):
    """pad peptide to target length with GSGG"""
    for seqrecord in SeqIO.parse(input, 'fasta'):
        padded = pad_ggsg(seqrecord.seq, length, terminus)
        pad_len = len(padded) - len(seqrecord)
        if pad_len > 0:
            output_title = '{}|{}-PADDED-{}'.format(seqrecord.id, terminus,
                                                    pad_len)
        else:
            output_title = seqrecord.id
        output_record = SeqRecord(padded, output_title, description='')
        SeqIO.write(output_record, output, 'fasta')


@cli.command()
@argument_input
@argument_output
@option('--prefix', '-p', help='DNA sequence to prefix each oligo')
def prefix(input, output, prefix):
    """add a prefix to each sequence"""
    for seqrecord in SeqIO.parse(input, 'fasta'):
        newseq = prefix + seqrecord.seq
        seqrecord.seq = newseq
        SeqIO.write(seqrecord, output, 'fasta')


@cli.command()
@argument_input
@argument_output
@option('--suffix', '-s', help='DNA sequence to suffix each oligo')
def suffix(input, output, suffix):
    """add a suffix to each sequence"""
    for seqrecord in SeqIO.parse(input, 'fasta'):
        newseq = seqrecord.seq + suffix
        seqrecord.seq = newseq
        SeqIO.write(seqrecord, output, 'fasta')


@cli.command()
@argument_input
@argument_output
@option('--codon-table', '-t', default='standard',
        help='ONLY STANDARD TABLE IMPLEMENTED')
@option('--codon-usage', '-u', default='ecoli', help='ONLY ECOLI IMPLEMENTED')
@option('--sampler', default='weighted', show_default=True,
        type=Choice(['weighted', 'uniform']), help='Codon sampling method')
@option('--codon-freq-threshold', type=float, default=None,
        help='Minimum codon frequency')
@option('--amber-only', is_flag=True, help='Use only amber stop codon')
def revtrans(input, output, codon_table, codon_usage, sampler,
             codon_freq_threshold, amber_only):
    """reverse translate amino acid sequences into DNA"""
    if sampler == 'weighted':
        usage = ecoli_codon_usage
        if codon_freq_threshold is not None:
            # TODO: this is hardcoded in and there's a leaky abstraction here
            table = standard_dna_table
            usage = zero_low_freq_codons(usage, table, codon_freq_threshold)
        if amber_only:
            usage = zero_non_amber_stops(usage)
        codon_sampler = FreqWeightedCodonSampler(usage=usage)
    elif sampler == 'uniform':
        codon_sampler = UniformCodonSampler()
    for seqrecord in SeqIO.parse(input, 'fasta'):
        dna_id = seqrecord.id
        dna_seq = reverse_translate(seqrecord.seq, codon_sampler)
        SeqIO.write(SeqRecord(dna_seq, dna_id, description=''), output, 'fasta')


@cli.command()
@argument_input
@argument_output
@option('--site', help='Site to remove (e.g., EcoRI, AGCCT); case sensitive')
@option('--clip-left', type=int, default=0,
        help='Number of bases to clip from start of sequence to get to CDS')
@option('--clip-right', type=int, default=0,
        help='Number of bases to clip from end of sequence to get to CDS')
@option('--codon-table', '-t', default='standard',
        help='ONLY STANDARD TABLE IMPLEMENTED')
@option('--codon-usage', '-u', default='ecoli', help='ONLY ECOLI IMPLEMENTED')
@option('--sampler', default='weighted', show_default=True,
        type=Choice(['weighted', 'uniform']), help='Codon sampling method')
@option('--codon-freq-threshold', type=float, default=None,
        help='Minimum codon frequency')
@option('--amber-only', is_flag=True, help='Use only amber stop codon')
def recodesite(input, output, site, clip_left, clip_right, codon_table,
               codon_usage, sampler, codon_freq_threshold, amber_only):
    """remove site from each sequence's CDS by recoding"""
    if sampler == 'weighted':
        usage = ecoli_codon_usage
        if codon_freq_threshold is not None:
            # TODO: this is hardcoded in and there's a leaky abstraction here
            table = standard_dna_table
            usage = zero_low_freq_codons(usage, table, codon_freq_threshold)
        if amber_only:
            usage = zero_non_amber_stops(usage)
        codon_sampler = FreqWeightedCodonSampler(usage=usage)
    elif sampler == 'uniform':
        codon_sampler = UniformCodonSampler()

    site = site2dna(site)
    # site is now a Bio.Seq.Seq

    for seqrecord in SeqIO.parse(input, 'fasta'):
        id_ = seqrecord.id
        cds_start = clip_left
        cds_end = len(seqrecord) - clip_right
        seq = recode_site_from_cds(seqrecord.seq, site, codon_sampler,
                                   cds_start, cds_end)
        SeqIO.write(SeqRecord(seq, id_, description=''), output, 'fasta')


@cli.command()
@argument_input
@argument_output
def x2ggsg(input, output):
    """replace stretches of Xs with Serine-Glycine linker (GGSG pattern)"""
    for seqrecord in SeqIO.parse(input, 'fasta'):
        replacement = x_to_ggsg(seqrecord.seq)
        if replacement != seqrecord.seq:
            output_title = '{}|{}'.format(seqrecord.id, 'withGSlinker')
        else:
            output_title = seqrecord.id
        output_record = SeqRecord(replacement, id=output_title, description='')
        SeqIO.write(output_record, output, 'fasta')


@cli.command()
@argument_input
@argument_output
def disambiguateaa(input, output):
    """replace ambiguous (IUPAC) AAs with unambiguous ones (e.g. Z => E/Q)

    B => DN, X => ACDEFGHIKLMNPQRSTVWY, Z => EQ, J => LI,
    U => C (selenocysteine), O => K (pyrrolysine)
    """
    for seqrecord in SeqIO.parse(input, 'fasta'):
        id_ = seqrecord.id
        ambig = seqrecord.seq
        for (i, unambig) in enumerate(disambiguate_iupac_aa(ambig)):
            if unambig != ambig:
                output_title = '{}|disambig_{}'.format(id_, i + 1)
            else:
                output_title = id_
            output_record = SeqRecord(unambig, output_title, description='')
            SeqIO.write(output_record, output, 'fasta')


@cli.command()
@argument_input
@option('--site', help='Site to find (e.g., EcoRI, AGCCT); case sensitive')
@option('--clip-left', type=int, default=0,
        help='Number of bases to clip from start of sequence to get to CDS')
@option('--clip-right', type=int, default=0,
        help='Number of bases to clip from end of sequence to get to CDS')
def findsite(input, site, clip_left, clip_right):
    """find locations of a site"""
    query = site2dna(site)
    for seqrecord in SeqIO.parse(input, 'fasta'):
        id_ = seqrecord.id
        start = clip_left
        end = len(seqrecord) - clip_right
        idx = seqrecord.seq[start:end].find(query)
        if idx >= 0:
            print('{}|{}|{}'.format(id_, site, idx + start), flush=True)


@cli.command()
@argument_input
def stats(input):
    """NOT IMPL'd: compute some sequence statistics"""
    pass
