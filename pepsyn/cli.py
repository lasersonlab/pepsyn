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

import sys
from math import ceil
from collections import Counter
from logging import captureWarnings

# biopython has a bunch of annoying warnings bc Seq comparisons changed
captureWarnings(True)

from click import (
    group, command, option, argument, File, Choice, Path, version_option)
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Data.CodonTable import standard_dna_table
from tqdm import tqdm, trange

from pepsyn import __version__
from pepsyn.operations import (
    reverse_translate, recode_site_from_cds, x_to_ggsg, disambiguate_iupac_aa,
    tile as tile_op, ctermpep as cterm_oligo, pad_ggsg)
from pepsyn.codons import (
    FreqWeightedCodonSampler, UniformCodonSampler, ecoli_codon_usage,
    zero_non_amber_stops, zero_low_freq_codons)
from pepsyn.util import site2dna


@group(context_settings={'help_option_names': ['-h', '--help']})
@version_option(__version__)
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
    for seqrecord in tqdm(SeqIO.parse(input, 'fasta'), desc='tile', unit='seq'):
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
    for seqrecord in tqdm(SeqIO.parse(input, 'fasta'), desc='ctermpep', unit='seq'):
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
    for seqrecord in tqdm(SeqIO.parse(input, 'fasta'), desc='revtrans', unit='seq'):
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
            tqdm.write('{}|{}|{}'.format(id_, site, idx + start))


@cli.command()
@argument_input
@argument_output
@option('-k', '--kmer-size', type=int, required=True, help='k-mer size')
@option('-t', '--tile-size', type=int, required=True, help='tile size')
@option('-c', '--kmer-cov', type=float, default=1.5,
        help='target avg k-mer coverage')
@option('-b', '--cterm-boost', type=float, default=1.,
        help='multiplier to increase representation of cterm tiles')
@option('-p', '--prefix', default='tile', help='tile size')
def greedykmercov(input, output, kmer_size, tile_size, kmer_cov, cterm_boost, prefix):
    """select protein tiles by maximizing k-mer coverage

    each tile is a fragment of an observed input ORF
    ORFs shorter than tile-size are still sampled
    ORFs shorter than kmer-size are always included in tile set
    """
    import networkx as nx
    from pepsyn.dbg import (
        gen_kmers, tiling_stats, orf_stats, seqrecords_to_dbg,
        sequence_incr_attr, sequence_setreduce_attr)

    # load data and build dbg
    orfs = {sr.id: sr for sr in SeqIO.parse(input, 'fasta')}
    with tqdm(desc='loading dbg') as pbar:
        dbg = seqrecords_to_dbg(
            orfs.values(), kmer_size, skip_short=True, tqdm=pbar)

    selected_tiles = set()

    # add ORFs shorter than kmer size to the selected tiles
    short_orf_names = {}
    for orf in orfs.values():
        if len(orf) < kmer_size:
            selected_tiles.add(str(orf.seq))
            short_orf_names[str(orf.seq)] = orf.id

    # process each graph component separately
    num_components = nx.number_weakly_connected_components(dbg)
    component_iter = nx.weakly_connected_component_subgraphs(dbg)
    for component in tqdm(component_iter, desc='components', total=num_components):
        # generate all candidate tiles
        cdss = set([cds for p in component.nodes_iter(True) for cds in p[1]['cds']])
        component_tiles = set()
        cterm_tiles = set()
        kmer_to_tile = {}
        for cds in cdss:
            cterm_tiles.add(str(orfs[cds].seq)[-tile_size:])
            for tile in gen_kmers(str(orfs[cds].seq), tile_size, yield_short=True):
                component_tiles.add(tile)
                for kmer in gen_kmers(tile, kmer_size):
                    kmer_to_tile.setdefault(kmer, set()).add(tile)

        def compute_tile_score(tile):
            weighted_multiplicity = 0
            for kmer in gen_kmers(tile, kmer_size, yield_short=True):
                if component.has_node(kmer) and component.node[kmer].get('weight', 0) == 0:
                    weighted_multiplicity += component.node[kmer]['multiplicity']
            return weighted_multiplicity

        # initialize tile scores
        tile_scores = {}
        for tile in component_tiles:
            tile_scores[tile] = compute_tile_score(tile)
            # give a boost to cterm tiles
            if tile in cterm_tiles:
                tile_scores[tile] += ceil(tile_size * cterm_boost)

        def update_tile_scores(tile):
            for kmer in gen_kmers(tile, kmer_size):
                for t in kmer_to_tile[kmer]:
                    tile_scores[t] = max(
                        0, tile_scores[t] - dbg.node[kmer]['multiplicity'])

        num_component_tiles = ceil(len(component) * kmer_cov / (tile_size - kmer_size + 1))
        for i in trange(num_component_tiles, desc='tile selection'):
            tile = max(component_tiles, key=tile_scores.get)
            if len(tile) < kmer_size:
                continue
            selected_tiles.add(tile)
            sequence_incr_attr(component, tile, kmer_size, 'weight')
            sequence_incr_attr(dbg, tile, kmer_size, 'weight')
            update_tile_scores(tile)

    for (k, v) in orf_stats(dbg, orfs.values(), tile_size):
        print(f'{k}\t{v}', file=sys.stderr)
    for (k, v) in tiling_stats(dbg, selected_tiles):
        print(f'{k}\t{v}', file=sys.stderr)

    # write out tiles and generate names
    for i, tile in enumerate(tqdm(selected_tiles, desc='writing', unit='tile')):
        if len(tile) < kmer_size:
            cdss = set([short_orf_names[tile]])
        else:
            cdss = sequence_setreduce_attr(dbg, tile, kmer_size, 'cds')
        cds = Counter(cdss).most_common(1)[0][0]
        print(f'>{prefix}{i:05d}|{cds}\n{tile}', file=output)


@cli.command()
@option('-p', '--tiles', type=Path(exists=True, dir_okay=False), required=True,
        help='input protein tiles fasta')
@option('-r', '--orfs', type=Path(exists=True, dir_okay=False), required=True,
        help='input ORFs fasta')
@option('-k', '--kmer-size', type=int, required=True, help='k-mer size')
@option('-t', '--tile-size', type=int, required=True, help='tile size')
def proteintilestats(tiles, orfs, kmer_size, tile_size):
    """compute some sequence statistics"""
    from pepsyn.dbg import (
        tiling_stats, orf_stats, fasta_to_dbg, sequence_incr_attr)

    with tqdm(desc='loading dbg') as pbar:
        dbg = fasta_to_dbg(orfs, kmer_size, skip_short=True, tqdm=pbar)
    orfs = list(SeqIO.parse(orfs, 'fasta'))

    # annotate DBG with tile weights
    tiles = [str(sr.seq) for sr in SeqIO.parse(tiles, 'fasta')]
    for tile in tiles:
        sequence_incr_attr(dbg, tile, kmer_size, 'weight')

    for (k, v) in orf_stats(dbg, orfs, tile_size):
        print(f'{k}\t{v}', file=sys.stderr)
    for (k, v) in tiling_stats(dbg, tiles):
        print(f'{k}\t{v}', file=sys.stderr)
