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

import os
import sys
from os.path import join as pjoin
from math import ceil, inf
from collections import Counter
from logging import captureWarnings

# biopython has a bunch of annoying warnings bc Seq comparisons changed
captureWarnings(True)

from click import (
    group, command, option, argument, File, Choice, Path, version_option,
    Abort, UsageError)
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Data.CodonTable import standard_dna_table
from tqdm import tqdm, trange
import yaml

from pepsyn import __version__
from pepsyn.operations import (
    reverse_translate, recode_site_from_cds, recode_sites_from_cds, x_to_ggsg,
    disambiguate_iupac_aa, tile as tile_op, ctermpep as cterm_oligo, pad_ggsg)
from pepsyn.codons import (
    FreqWeightedCodonSampler, UniformCodonSampler, ecoli_codon_usage,
    zero_non_amber_stops, zero_low_freq_codons)
from pepsyn.util import site2dna, sliding_window, readfq


def print_fasta(sr, out):
    print('>{}\n{}'.format(sr.id, str(sr.seq)), file=out)


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
    for (name, seq, qual) in tqdm(readfq(input), desc='tile', unit='seq'):
        for (start, end, t) in tile_op(seq, length, overlap):
                output_title = f'{name}|{start}-{end}'
                print(f'>{output_title}\n{t}', file=output)


@cli.command()
@argument_input
@argument_output
@option('--length', '-l', type=int, help='Length of output C-terminal oligo')
@option('--add-stop', '-s', is_flag=True, help='Add a stop codon to peptide')
def ctermpep(input, output, length, add_stop):
    """extract C-terminal peptide (AA alphabets)

    will return entirety of seqs shorter than length
    """
    for (name, seq, qual) in tqdm(readfq(input), desc='ctermpep', unit='seq'):
        oligo = cterm_oligo(seq, length, add_stop=add_stop)
        output_title = f'{name}|CTERM'
        if add_stop:
            output_title = f'{output_title}|STOP'
        print(f'>{output_title}\n{oligo}', file=output)


@cli.command()
@argument_input
@argument_output
def stripstop(input, output):
    """strip stop codons from end of protein sequence"""
    for (name, seq, qual) in readfq(input):
        seq = seq.rstrip('*')
        print(f'>{name}\n{seq}', file=output)


@cli.command()
@argument_input
@argument_output
def filterstop(input, output):
    """filter out sequences that contain stop codons (*)"""
    for (name, seq, qual) in readfq(input):
        if '*' in seq:
            continue
        seq = seq.rstrip('*')
        print(f'>{name}\n{seq}', file=output)


@cli.command()
@argument_input
@argument_output
@option('-m', '--min-len', type=int, help='Min length of sequence to keep')
@option('-M', '--max-len', type=int, help='Max length of sequence to keep')
def filterlen(input, output, min_len, max_len):
    """filter sequences of a given length"""
    if min_len is None:
        min_len = 0
    if max_len is None:
        max_len = inf
    for (name, seq, qual) in readfq(input):
        if len(seq) < min_len:
            continue
        if len(seq) > max_len:
            continue
        print(f'>{name}\n{seq}', file=output)


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
    for (name, seq, qual) in readfq(input):
        padded = pad_ggsg(seq, length, terminus)
        pad_len = len(padded) - len(seq)
        if pad_len > 0:
            output_title = f'{name}|{terminus}-PADDED-{pad_len}'
        else:
            output_title = name
        print(f'>{output_title}\n{padded}', file=output)


@cli.command()
@argument_input
@argument_output
@option('--prefix', '-p', help='DNA sequence to prefix each oligo')
def prefix(input, output, prefix):
    """add a prefix to each sequence"""
    for (name, seq, qual) in readfq(input):
        newseq = prefix + seq
        print(f'>{name}\n{newseq}', file=output)


@cli.command()
@argument_input
@argument_output
@option('--suffix', '-s', help='DNA sequence to suffix each oligo')
def suffix(input, output, suffix):
    """add a suffix to each sequence"""
    for (name, seq, qual) in readfq(input):
        newseq = seq + suffix
        print(f'>{name}\n{newseq}', file=output)


@cli.command()
@argument_input
@argument_output
@option('--left', '-l', default=0, help='number of bases to clip from left')
@option('--right', '-r', default=0, help='number of bases to clip from right')
def clip(input, output, left, right):
    """clip bases from the ends of each sequence"""
    for (name, seq, qual) in readfq(input):
        stop = len(seq) - right
        print(f'>{name}\n{seq[left:stop]}', file=output)


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
        print_fasta(SeqRecord(dna_seq, dna_id, description=''), output)


@cli.command()
@argument_input
@argument_output
@option('--site', multiple=True,
        help='Site to remove (e.g., EcoRI, AGCCT); case sens; allow multiple')
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

    sites = [site2dna(s) for s in site]
    # sites is now a list[Bio.Seq.Seq]

    for seqrecord in SeqIO.parse(input, 'fasta'):
        id_ = seqrecord.id
        cds_start = clip_left
        cds_end = len(seqrecord) - clip_right
        seq = recode_sites_from_cds(seqrecord.seq, sites, codon_sampler,
                                    cds_start, cds_end)
        print_fasta(SeqRecord(seq, id_, description=''), output)


@cli.command()
@argument_input
@argument_output
def x2ggsg(input, output):
    """replace stretches of Xs with Serine-Glycine linker (GGSG pattern)"""
    for (name, seq, qual) in readfq(input):
        replacement = x_to_ggsg(seq)
        if replacement != seq:
            output_title = f'{name}|withGSlinker'
        else:
            output_title = name
        print(f'>{output_title}\n{replacement}', file=output)


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
            print_fasta(output_record, output)


@cli.command()
@argument_input
@option('--site', help='Site to find (e.g., EcoRI, AGCCT); case sensitive')
@option('--clip-left', type=int, default=0,
        help='Number of bases to clip from start of sequence to get to CDS')
@option('--clip-right', type=int, default=0,
        help='Number of bases to clip from end of sequence to get to CDS')
def findsite(input, site, clip_left, clip_right):
    """find locations of a site"""
    query = str(site2dna(site))
    for (name, seq, qual) in readfq(input):
        start = clip_left
        end = len(seq) - clip_right
        idx = seq[start:end].find(query)
        if idx >= 0:
            print(f'{name}|{site}|{idx + start}', file=output)


@cli.command()
@argument_input
@argument('output_path', type=Path(exists=False))
@option('-k', '--kmer-size', type=int, required=True, help='k-mer size')
def builddbg(input, output_path, kmer_size):
    """build a de bruijn graph

    ignores input sequences shorter than the kmer size

    output path with .gz will write compressed data
    """
    try:
        import networkx as nx
        from pepsyn.dbg import fasta_handle_to_dbg
    except ImportError:
        raise Abort('builddbg requires NetworkX')

    with tqdm(desc='building dbg') as pbar:
        dbg = fasta_handle_to_dbg(input, kmer_size, tqdm=pbar,
                                  ignore_short=True)
    nx.write_gpickle(dbg, output_path)


@cli.command()
@argument_input
@argument_output
@option('-t', '--tile-size', type=int, required=True, help='tile size')
@option('-d', '--dbg-path', type=Path(exists=True, dir_okay=False),
        required=True,
        help='path to prebuilt de bruijn graph (pickle file; .gz okay)')
@option('-c', '--kmer-cov', type=float,
        help='target avg k-mer coverage (incompatible with -n)')
@option('-n', '--num-tiles', type=int,
        help='number of output tiles (incompatible with -c)')
@option('-p', '--preselected-tiles-path', type=File('r'),
        help='fasta file containing previously selected tiles')
def greedykmercov(input, output, tile_size, dbg_path, kmer_cov, num_tiles,
                  preselected_tiles_path):
    """select protein tiles by maximizing k-mer coverage on de bruijn graph

    each tile is a fragment of an observed input ORF
    ORFS shorter than tile-size are sampled, but
    ORFs shorter than kmer-size are ignored
    """
    # test input/context
    try:
        import networkx as nx
        from pepsyn.dbg import gen_kmers, setreduce_attr, sum_attr
    except ImportError:
        raise Abort('greedykmercov requires NetworkX')
    try:
        import numpy as np
    except ImportError:
        raise Abort('greedykmercov requires NumPy')
    if kmer_cov and num_tiles:
        raise UsageError('Set -c/--kmer-cov OR -n/--num-tiles but not both')
    if not kmer_cov and not num_tiles:
        raise UsageError('Must set one of -c/--kmer-cov OR -n/--num-tiles')

    # load orfs
    orfs = {name: seq for (name, seq, qual) in readfq(input)}

    # load dbg
    dbg = nx.read_gpickle(dbg_path)
    kmer_size = len(next(iter(dbg)))
    if kmer_size > tile_size:
        raise UsageError('kmer-size > tile_size')
    kmers_remaining = len(dbg)
    num_components = nx.number_weakly_connected_components(dbg)
    if num_tiles:
        tiles_remaining = num_tiles

    # load preselected tiles
    preselected_tiles = [seq for (name, seq, qual) in readfq(preselected_tiles_path)]
    preselected_kmer_counts = Counter(
        [kmer for tile in preselected_tiles for kmer in gen_kmers(tile, kmer_size, yield_short=True)])

    # process each graph component separately
    component_iter = tqdm(nx.weakly_connected_components(dbg), unit='comp',
                          desc='dbg components', total=num_components)
    for component in component_iter:
        component_orfs = setreduce_attr(dbg, component, 'orf')

        # generate all candidate tiles
        tile_to_name = {}
        for name in tqdm(component_orfs, desc='generating tiles'):
            # special case short orfs
            if len(orfs[name]) < tile_size:
                tile_to_name.setdefault(orfs[name], []).append((name, 0, len(orfs[name])))
            for (i, j, tile) in tile_op(orfs[name], tile_size, tile_size - 1):
                tile_to_name.setdefault(tile, []).append((name, i, j))
        candidate_tiles = list(tile_to_name.keys())

        # generate init tile scores
        tile_scores = []
        tile_lens = []
        kmer_to_idxs = {}
        for idx, tile in enumerate(tqdm(candidate_tiles, desc='init tile scores')):
            score = 0
            for kmer in set(gen_kmers(tile, kmer_size)):
                score += dbg.nodes[kmer]['multiplicity']
                kmer_to_idxs.setdefault(kmer, set()).add(idx)
            tile_scores.append(score / len(tile))
            tile_lens.append(len(tile))
        tile_scores = np.ma.asarray(tile_scores)
        tile_scores.harden_mask()
        tile_lens = np.asarray(tile_lens)

        # update tile scores with previously selected tiles
        for kmer in set(preselected_kmer_counts.keys()) & set(kmer_to_idxs.keys()):
            idxs = list(kmer_to_idxs[kmer])
            tile_scores.data[idxs] -= (
                (preselected_kmer_counts[kmer] * dbg.nodes[kmer]['multiplicity']) / len(tile))

        # set number of tiles for this component
        if kmer_cov:
            num_component_tiles = ceil(len(component) * kmer_cov / (tile_size - kmer_size + 1))
        if num_tiles:
            num_component_tiles = ceil(len(component) / kmers_remaining * tiles_remaining)
            kmers_remaining -= len(component)
            tiles_remaining -= num_component_tiles

        # choose tiles
        for _ in trange(num_component_tiles, desc='choosing tiles'):
            idx = tile_scores.argmax()
            tile_scores[idx] = np.ma.masked
            tile = candidate_tiles[idx]

            # write tile
            name, i, j = tile_to_name[tile][0]
            nterm = '|NTERM' if dbg.nodes[tile[:kmer_size]].get('start_node', False) else ''
            cterm = '|CTERM' if dbg.nodes[tile[-kmer_size:]].get('end_node', False) else ''
            print(f'>{name}|{i}-{j}{nterm}{cterm}\n{tile}', file=output)

            # update tile scores
            for kmer in set(gen_kmers(tile, kmer_size)):
                idxs = list(kmer_to_idxs[kmer])
                tile_scores.data[idxs] -= dbg.nodes[kmer]['multiplicity'] / tile_lens[idxs]


@cli.command()
@option('-d', '--dbg-path', type=Path(exists=True, dir_okay=False), required=True,
        help='pickled dbg')
@option('-p', '--tiles', type=Path(exists=True, dir_okay=False), required=True,
        help='input protein tiles fasta')
@option('-r', '--orfs', type=Path(exists=True, dir_okay=False), required=True,
        help='input ORFs fasta')
@option('-o', '--output', type=File('w'), default=sys.stdout,
        help='output YAML path (default stdout)')
def dbgtilesummary(dbg_path, tiles, orfs, output):
    """compute some sequence statistics"""
    try:
        import networkx as nx
        from pepsyn.dbg import dbg_orf_tile_stats
    except ImportError:
        raise Abort('dbgtilesummary requires NetworkX')

    dbg = nx.read_gpickle(dbg_path)
    with open(orfs, 'r') as ip:
        orfs = [seq for (name, seq, qual) in readfq(ip)]
    with open(tiles, 'r') as ip:
        tiles = [seq for (name, seq, qual) in readfq(ip)]

    stats = dbg_orf_tile_stats(dbg, orfs, tiles)

    print(yaml.dump(stats), file=output)
