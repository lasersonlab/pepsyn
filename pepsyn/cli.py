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
    """extract C-terminal peptide (AA alphabets)"""
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
@argument_output
@option('-k', '--kmer-size', type=int, required=True, help='k-mer size')
@option('-t', '--tile-size', type=int, required=True, help='tile size')
@option('-c', '--kmer-cov', type=float,
        help='target avg k-mer coverage (incompatible with -n)')
@option('-n', '--num-tiles', type=int,
        help='number of output tiles (incompatible with -c)')
@option('-d', '--nterm-boost', type=float, default=1.,
        help='multiplier to increase representation of nterm tiles')
@option('-b', '--cterm-boost', type=float, default=1.,
        help='multiplier to increase representation of cterm tiles')
@option('-p', '--prefix', default='tile', help='tile size')
@option('-u', '--unweighted', is_flag=True,
        help='use unweighted k-mer coverage')
def greedykmercov(input, output, kmer_size, tile_size, kmer_cov, num_tiles,
                  nterm_boost, cterm_boost, prefix, unweighted):
    """select protein tiles by maximizing k-mer coverage

    each tile is a fragment of an observed input ORF
    ORFs shorter than tile-size are still sampled
    ORFs shorter than kmer-size are always included in tile set
    """
    import networkx as nx
    import numpy as np
    from pepsyn.dbg import (
        gen_kmers, tiling_stats, orf_stats, seqrecords_to_dbg, setreduce_attr,
        seq_to_path, path_to_seq, incr_attr)

    if kmer_cov and num_tiles:
        print('Set -c/--kmer-cov OR -n/--num-tiles but not both',
            file=sys.stderr)
        sys.exit(1)
    if not kmer_cov and not num_tiles:
        print('Must set one of -c/--kmer-cov OR -n/--num-tiles',
            file=sys.stderr)
        sys.exit(1)

    # load data and build dbg
    orfs = {sr.id: sr for sr in SeqIO.parse(input, 'fasta')}
    with tqdm(desc='building dbg') as pbar:
        dbg = seqrecords_to_dbg(
            orfs.values(), kmer_size, skip_short=True, tqdm=pbar)

    selected_tiles = set()

    # add ORFs shorter than kmer size to the selected tiles
    short_orf_names = {}
    for orf in orfs.values():
        if len(orf) < kmer_size:
            tile = str(orf.seq)
            selected_tiles.add(tile)
            short_orf_names[tile] = orf.id

    # process each graph component separately
    num_kmers = len(dbg)
    num_components = nx.number_weakly_connected_components(dbg)
    component_iter = nx.weakly_connected_components(dbg)
    for component in tqdm(component_iter, desc='tile selection', total=num_components):
        component_cdss = setreduce_attr(dbg, component, 'cds')

        # generate all candidate tiles/paths
        component_paths = []
        nterm_paths = []
        cterm_paths = []
        for cds in tqdm(component_cdss, desc='generating tiles'):
            orf_path = seq_to_path(str(orfs[cds].seq), kmer_size)

            if len(orfs[cds]) < tile_size:
                paths = [orf_path]
            else:
                paths = list(sliding_window(tile_size - kmer_size + 1, orf_path))

            nterm_paths.append(paths[0])
            cterm_paths.append(paths[-1])
            for path in paths:
                component_paths.append(path)
        component_paths = list(set(component_paths))
        nterm_paths = frozenset(nterm_paths)
        cterm_paths = frozenset(cterm_paths)

        kmer_to_idxs = {}
        for i, path in enumerate(tqdm(component_paths, desc='kmer mapping')):
            for kmer in path:
                kmer_to_idxs.setdefault(kmer, []).append(i)

        # initialize path scores
        path_scores = []
        for path in tqdm(component_paths, desc='init scores'):
            if unweighted:
                score = len(path)
            else:
                score = sum(dbg.node[kmer]['multiplicity'] for kmer in path)
            # give a boost to cterm or nterm tiles
            if path in nterm_paths:
                score += ceil(tile_size * nterm_boost)
            if path in cterm_paths:
                score += ceil(tile_size * cterm_boost)
            path_scores.append(score)
        path_scores = np.ma.asarray(path_scores)
        path_scores.harden_mask()

        # set number of tiles for this component
        if kmer_cov:
            num_component_tiles = ceil(len(component) * kmer_cov / (tile_size - kmer_size + 1))
        if num_tiles:
            num_component_tiles = ceil(len(component) / num_kmers * num_tiles)

        # choose tiles
        for _ in trange(num_component_tiles, desc='choosing tiles'):
            i = path_scores.argmax()
            path_scores[i] = np.ma.masked
            path = component_paths[i]
            selected_tiles.add(path_to_seq(path))
            for kmer in path:
                # update weight (coverage) on graph
                incr_attr(dbg, kmer, 'weight')
                # update path scores
                idxs = list(set(kmer_to_idxs[kmer]))
                if unweighted:
                    path_scores.data[idxs] = path_scores.data[idxs] - 1
                else:
                    path_scores.data[idxs] = path_scores.data[idxs] - dbg.node[kmer]['multiplicity']
                zeros = path_scores.data < 0
                path_scores.data[zeros] = 0

    for (k, v) in orf_stats(dbg, orfs.values(), tile_size):
        print(f'{k}\t{v}', file=sys.stderr)
    for (k, v) in tiling_stats(dbg, selected_tiles):
        print(f'{k}\t{v}', file=sys.stderr)

    # write out tiles and generate names
    for i, tile in enumerate(selected_tiles):
        nterm = ''
        cterm = ''
        if len(tile) < kmer_size:
            cdss = set([short_orf_names[tile]])
            nterm = '|NTERM'
            cterm = '|CTERM'
        else:
            cdss = setreduce_attr(dbg, gen_kmers(tile, kmer_size), 'cds')
            if dbg.node[tile[:kmer_size]].get('nterm', False):
                nterm = '|NTERM'
            if dbg.node[tile[-kmer_size:]].get('end_node', False):
                cterm = '|CTERM'
        cds = Counter(cdss).most_common(1)[0][0]
        print(f'>{prefix}{i:05d}|{cds}{nterm}{cterm}\n{tile}', file=output)


@cli.command()
@option('-p', '--tiles', type=Path(exists=True, dir_okay=False), required=True,
        help='input protein tiles fasta')
@option('-r', '--orfs', type=Path(exists=True, dir_okay=False), required=True,
        help='input ORFs fasta')
@option('-k', '--kmer-size', type=int, required=True, help='k-mer size')
@option('-t', '--tile-size', type=int, required=True, help='tile size')
def proteintilestats(tiles, orfs, kmer_size, tile_size):
    """compute some sequence statistics

    skips tiles shorter than kmer-size
    """
    from pepsyn.dbg import (
        tiling_stats, orf_stats, fasta_to_dbg, sequence_incr_attr)

    with tqdm(desc='building dbg') as pbar:
        dbg = fasta_to_dbg(orfs, kmer_size, skip_short=True, tqdm=pbar)
    orfs = list(SeqIO.parse(orfs, 'fasta'))

    # annotate DBG with tile weights
    tiles = [str(sr.seq) for sr in SeqIO.parse(tiles, 'fasta')]
    for tile in tiles:
        if len(tile) < kmer_size:
            continue
        sequence_incr_attr(dbg, tile, kmer_size, 'weight')

    for (k, v) in orf_stats(dbg, orfs, tile_size):
        print(f'{k}\t{v}', file=sys.stderr)
    for (k, v) in tiling_stats(dbg, tiles):
        print(f'{k}\t{v}', file=sys.stderr)
