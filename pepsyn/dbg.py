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
    dbg.node[kmers[0]]['nterm'] = True
    dbg.node[kmers[-1]]['cterm'] = True


def orf_stats(dbg, orfs, tile_size):
    stats = []
    kmer_size = len(dbg.nodes()[0])
    stats.append(('kmer size', kmer_size))
    stats.append(('num ORFs', len(orfs)))
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
    stats.append(('num kmers covered', graph_num_pos_attr(dbg, 'weight')))
    stats.append(
        ('frac kmers covered', graph_num_pos_attr(dbg, 'weight') / len(dbg)))
    stats.append(
        ('avg kmer coverage', graph_sum_attr(dbg, 'weight') / len(dbg)))
    stats.append(('max kmer coverage', graph_max_attr(dbg, 'weight')))
    weighted_cov = sum([d['multiplicity'] for d in dbg.node.values() if d.get('weight', 0) > 0]) / graph_sum_attr(dbg, 'multiplicity')
    stats.append(('multiplicity-weighted kmer coverage', weighted_cov))
    cterm_cov = sum([1 for d in dbg.node.values() if (d.get('weight', 0) and d.get('cterm', False))]) / graph_sum_attr(dbg, 'cterm')
    stats.append(('c-term kmer coverage', cterm_cov))
    return stats




import graph_tool.all as gt

def gt_build_dbg(fasta_file):
    k = 15
    kmer_to_index = {}
    dbg = gt.Graph()
    dbg.vp.kmer = dbg.new_vertex_property('string')
    dbg.vp.nterm = dbg.new_vertex_property('bool')
    dbg.vp.cterm = dbg.new_vertex_property('bool')
    dbg.vp.multiplicity = dbg.new_vertex_property('int32_t')
    dbg.vp.cds = dbg.new_vertex_property('vector<string>')
    for sr in tqdm(SeqIO.parse(fasta_file, 'fasta')):
        if len(sr) < k:
            continue
        kmer_path = list(gen_kmers(str(sr.seq), k))
        for kmer in kmer_path:
            if kmer not in kmer_to_index:
                v = dbg.add_vertex()
                kmer_to_index[kmer] = int(v)
                dbg.vp.kmer[v] = kmer
                dbg.vp.multiplicity[v] = 0
                dbg.vp.cds[v] = []
                dbg.vp.nterm[v] = False
                dbg.vp.cterm[v] = False
            else:
                v = dbg.vertex(kmer_to_index[kmer])
            dbg.vp.multiplicity[v] += 1
            dbg.vp.cds[v].append(sr.id)
        dbg.vp.nterm[kmer_to_index[kmer_path[0]]] = True
        dbg.vp.cterm[kmer_to_index[kmer_path[-1]]] = True

        for (kmer1, kmer2) in sliding_window(2, kmer_path):
            dbg.add_edge(kmer_to_index[kmer1], kmer_to_index[kmer2])
    return dbg
