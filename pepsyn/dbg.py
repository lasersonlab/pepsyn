# Copyright 2017 Uri Laserson
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


from math import ceil

import networkx as nx
from Bio import SeqIO


# some utilities for working with graphs


def incr_attr(graph, node, attr, amt=1):
    graph.node[node][attr] = graph.node[node].get(attr, 0) + amt


def zero_attr(graph, node, attr):
    graph.node[node][attr] = 0


def graph_zero_attr(graph, attr):
    for node in graph.nodes():
        graph.node[node][attr] = 0


def sum_attr(graph, nodes, attr):
    return sum([graph.node[node].get(attr, 0) for node in nodes])


def setreduce_attr(graph, nodes, attr):
    result = set()
    for node in nodes:
        result |= set(graph.node.get(node, {}).get(attr, []))
    return result


def max_attr(graph, nodes, attr):
    return max([graph.node[node].get(attr, 0) for node in nodes])


def graph_sum_attr(graph, attr):
    return sum([d.get(attr, 0) for d in graph.node.values()])


def graph_max_attr(graph, attr):
    return max([d.get(attr, 0) for d in graph.node.values()])


def graph_num_pos_attr(graph, attr):
    num_pos = 0
    for d in graph.node.values():
        if d.get(attr, 0) > 0:
            num_pos += 1
    return num_pos


# some utilities for working with sequences on a de bruijn graph


def seq_to_path(s, k):
    return tuple(gen_kmers(s, k))


def path_to_seq(p):
    return ''.join([p[0]] + [n[-1] for n in p[1:]])


def sequence_sum_attr(dbg, s, k, attr):
    return sum_attr(dbg, gen_kmers(s, k), attr)


def sequence_incr_attr(dbg, s, k, attr, amt=1):
    for kmer in gen_kmers(s, k):
        if dbg.has_node(kmer):
            incr_attr(dbg, kmer, attr, amt)


def sequence_setreduce_attr(dbg, s, k, attr):
    result = set()
    for kmer in gen_kmers(s, k):
        result |= set(dbg.node.get(kmer, {}).get(attr, []))
    return result


def gen_kmers(s, k, yield_short=False):
    if k <= 0:
        raise ValueError('k-mer len must be positive')
    if len(s) < k:
        if yield_short:
            yield s
        else:
            raise ValueError('k-mer len longer than string')
    for i in range(len(s) - k + 1):
        yield s[i:i + k]


def seqrecords_to_dbg(seqrecords, k, tqdm=None, skip_short=False):
    dbg = nx.DiGraph()
    for sr in seqrecords:
        if tqdm is not None:
            tqdm.update()
        if skip_short and len(sr) < k:
            continue
        update_debruijn_graph(dbg, sr, k)
    return dbg


def fasta_to_dbg(fasta_file, k, tqdm=None, skip_short=False):
    return seqrecords_to_dbg(
        SeqIO.parse(fasta_file, 'fasta'), k, tqdm=tqdm, skip_short=skip_short)


def update_debruijn_graph(dbg, sr, k):
    """update de bruijn graph with a new sequence

    dbg is nx.DiGraph to mutate
    sr is Bio.SequenceRecord
    k is int k-mer size
    """
    s = str(sr.seq)
    if len(s) < k:
        raise ValueError('len(record) < k:\n{}'.format(str(sr)))
    kmers = list(gen_kmers(s, k))
    # define the graph structure
    if len(kmers) == 1:
        dbg.add_nodes_from(kmers)
    else:
        dbg.add_path(kmers)
    # add cds labels to nodes and count "multiplicity"
    for kmer in kmers:
        dbg.node[kmer].setdefault('cds', []).append(sr.id)
        dbg.node[kmer]['multiplicity'] = dbg.node[kmer].get('multiplicity', 0) + 1
    # annotate N/C-terminal nodes
    dbg.node[kmers[0]]['start_node'] = True
    dbg.node[kmers[-1]]['end_node'] = True


def orf_stats(dbg, orfs, tile_size):
    stats = []
    kmer_size = len(dbg.nodes()[0])
    stats.append(('kmer size', kmer_size))
    stats.append(('num ORFs', len(orfs)))
    stats.append(('total ORF residues', sum([len(orf) for orf in orfs])))
    stats.append(
        ('num ORFs smaller than tile size',
         len(list(filter(lambda x: len(x) < tile_size, orfs)))))
    stats.append(
        ('num ORFs smaller than k-mer size',
         len(list(filter(lambda x: len(x) < kmer_size, orfs)))))
    stats.append(('total ORF-observed kmers', len(dbg)))
    stats.append(('max theoretical kmers per tile', tile_size - kmer_size + 1))
    stats.append(
        ('min theoretical tiles for perfect kmer cov',
         ceil(len(dbg) / (tile_size - kmer_size + 1))))
    stats.append(
        ('approx num tiles in naive 1x tiling',
         sum([ceil(len(orf) / tile_size) for orf in orfs])))
    multiplicities = [d['multiplicity'] for d in dbg.node.values()]
    stats.append(('num multiplicity-1 kmers', multiplicities.count(1)))
    stats.append(
        ('avg kmer multiplicity', sum(multiplicities) / len(multiplicities)))
    stats.append(('max kmer multiplicity', max(multiplicities)))
    stats.append(
        ('num components', nx.number_weakly_connected_components(dbg)))
    return stats


def tiling_stats(dbg, tiles):
    stats = []
    stats.append(('tile size', len(next(iter(tiles)))))
    stats.append(('num tiles', len(tiles)))
    stats.append(('total tiling residues', sum([len(tile) for tile in tiles])))
    stats.append(('num kmers covered', graph_num_pos_attr(dbg, 'weight')))
    stats.append(
        ('frac kmers covered', graph_num_pos_attr(dbg, 'weight') / len(dbg)))
    stats.append(
        ('avg kmer coverage', graph_sum_attr(dbg, 'weight') / len(dbg)))
    stats.append(('max kmer coverage', graph_max_attr(dbg, 'weight')))
    weighted_cov = sum([d['multiplicity'] for d in dbg.node.values() if d.get('weight', 0) > 0]) / graph_sum_attr(dbg, 'multiplicity')
    stats.append(('multiplicity-weighted kmer coverage', weighted_cov))
    nterm_cov = sum([1 for d in dbg.node.values() if (d.get('weight', 0) and d.get('start_node', False))]) / graph_sum_attr(dbg, 'start_node')
    stats.append(('n-term kmer coverage', nterm_cov))
    cterm_cov = sum([1 for d in dbg.node.values() if (d.get('weight', 0) and d.get('end_node', False))]) / graph_sum_attr(dbg, 'end_node')
    stats.append(('c-term kmer coverage', cterm_cov))
    return stats
