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
import numpy as np
from tqdm import tqdm

from pepsyn.util import compute_int_hist, readfq

# some utilities for working with graphs


def incr_attr(graph, node, attr, amt=1):
    graph.nodes[node][attr] = graph.nodes[node].get(attr, 0) + amt


def zero_attr(graph, node, attr):
    graph.nodes[node][attr] = 0


def graph_zero_attr(graph, attr):
    for node in graph:
        graph.nodes[node][attr] = 0


def sum_attr(graph, nodes, attr):
    return sum([graph.nodes[node].get(attr, 0) for node in nodes])


def setreduce_attr(graph, nodes, attr):
    result = set()
    for node in nodes:
        result |= set(graph.nodes.get(node, {}).get(attr, set()))
    return result


def max_attr(graph, nodes, attr):
    return max([graph.nodes[node].get(attr, 0) for node in nodes])


def graph_sum_attr(graph, attr):
    return sum([v for (_, v) in graph.nodes(data=attr, default=0)])


def graph_max_attr(graph, attr):
    return max([v for (_, v) in graph.nodes(data=attr, default=0)])


def graph_num_pos_attr(graph, attr):
    num_pos = 0
    for (_, v) in graph.nodes(data=attr, default=0):
        if v > 0:
            num_pos += 1
    return num_pos


# some utilities for working with sequences on a de bruijn graph


def seq_to_path(s, k):
    return tuple(gen_kmers(s, k))


def path_to_seq(p):
    return "".join([p[0]] + [n[-1] for n in p[1:]])


def sequence_sum_attr(dbg, s, k, attr):
    return sum_attr(dbg, gen_kmers(s, k), attr)


def sequence_incr_attr(dbg, s, k, attr, amt=1):
    for kmer in gen_kmers(s, k):
        if dbg.has_node(kmer):
            incr_attr(dbg, kmer, attr, amt)


def sequence_setreduce_attr(dbg, s, k, attr):
    result = set()
    for kmer in gen_kmers(s, k):
        result |= set(dbg.nodes.get(kmer, {}).get(attr, []))
    return result


def gen_kmers(s, k, yield_short=False):
    if k <= 0:
        raise ValueError("k-mer len must be positive")
    if len(s) < k:
        if yield_short:
            yield s
        else:
            raise ValueError("k-mer len longer than string")
    for i in range(len(s) - k + 1):
        yield s[i : i + k]


def update_debruijn_graph(dbg, k, seq, name=None):
    """update de bruijn graph with a new sequence

    dbg is nx.DiGraph to mutate
    k is int k-mer size
    seq is str
    name is optional str
    """
    if len(seq) < k:
        raise ValueError(f"len(seq) < k:\n{name}\n{seq}")
    # define the graph structure
    kmers = list(gen_kmers(seq, k))
    if len(kmers) == 1:
        dbg.add_nodes_from(kmers)
    else:
        dbg.add_path(kmers)
    # add orf labels to nodes and count "multiplicity"
    for kmer in kmers:
        if name:
            dbg.nodes[kmer].setdefault("orf", set()).add(name)
        dbg.nodes[kmer]["multiplicity"] = dbg.nodes[kmer].get("multiplicity", 0) + 1
    # annotate N/C-terminal nodes
    dbg.nodes[kmers[0]]["start_node"] = True
    dbg.nodes[kmers[-1]]["end_node"] = True


def name_seq_pairs_to_dbg(pairs, k, tqdm=None, ignore_short=False):
    """build de bruijn graph from sequences

    pairs is an iter of (name, seq) pairs (iter[(str, str)]); name can be None
    """
    dbg = nx.DiGraph()
    for (name, seq) in pairs:
        if tqdm is not None:
            tqdm.update()
        if ignore_short and len(seq) < k:
            continue
        update_debruijn_graph(dbg, k, seq, name)
    return dbg


def seqrecords_to_dbg(seqrecords, k, tqdm=None, ignore_short=False):
    return name_seq_pairs_to_dbg(
        ((sr.id, str(sr.seq)) for sr in seqrecords), k, tqdm, ignore_short
    )


def fasta_handle_to_dbg(fasta_handle, k, tqdm=None, ignore_short=False):
    return name_seq_pairs_to_dbg(
        ((name, seq) for (name, seq, qual) in readfq(fasta_handle)),
        k,
        tqdm,
        ignore_short,
    )


def fasta_file_to_dbg(fasta_file, k, tqdm=None, ignore_short=False):
    with open(fasta_file, "r") as ip:
        return fasta_handle_to_dbg(ip, k, tqdm, ignore_short)


def dbg_stats(dbg, orfs, tiles, cov_attr="coverage"):
    """compute de bruijn graph stats

    dbg is a de bruijn graph
    orfs and tiles are lists of seqs (strings)

    NOTE: this operation mutates the in-memory dbg by computing coverage of
    each kmer into the cov_attr attribute

    some kmers in the tiles may not be found in the dbg
    """
    tile_lens = np.asarray([len(t) for t in tiles])
    orf_lens = np.asarray([len(o) for o in orfs])
    kmer_size = len(next(iter(dbg)))
    tile_size = int(round(np.median(tile_lens)).tolist())

    # compute coverage of tiles on dbg
    for tile in tqdm(tiles, desc="computing coverage"):
        for kmer in gen_kmers(tile, kmer_size, yield_short=True):
            if dbg.has_node(kmer):
                dbg.nodes[kmer][cov_attr] = dbg.nodes[kmer].get(cov_attr, 0) + 1

    coverages = np.asarray([cov for (_, cov) in dbg.nodes(data=cov_attr, default=0)])
    multiplicities = np.asarray([mult for (_, mult) in dbg.nodes(data="multiplicity")])
    nterms = np.asarray([s if s else False for (_, s) in dbg.nodes(data="start_node")])
    cterms = np.asarray([e if e else False for (_, e) in dbg.nodes(data="end_node")])

    assert len(coverages) == len(dbg)
    assert len(multiplicities) == len(dbg)
    assert graph_sum_attr(dbg, "multiplicity") == multiplicities.sum()

    stats = {}
    stats["kmer_size"] = kmer_size
    stats["num_orfs_smaller_than_kmer_size"] = (orf_lens < kmer_size).sum().tolist()
    stats["num_observed_kmers"] = len(dbg)
    stats["max_theor_kmers_per_tile"] = tile_size - kmer_size + 1
    stats["min_theor_tiles_1x_kmer_cov"] = ceil(len(dbg) / (tile_size - kmer_size + 1))
    stats["num_multiplicity_1_kmers"] = (multiplicities == 1).sum().tolist()
    stats["num_multiplicity_gt1_kmers"] = (multiplicities > 1).sum().tolist()
    stats["avg_kmer_multiplicity"] = multiplicities.mean().tolist()
    stats["max_kmer_multiplicity"] = multiplicities.max().tolist()
    stats["num_dbg_components"] = nx.number_weakly_connected_components(dbg)
    stats["num_kmers_covered"] = (coverages > 0).sum().tolist()
    stats["num_kmers_missed"] = (coverages == 0).sum().tolist()
    stats["frac_kmers_covered"] = (coverages > 0).sum().tolist() / len(coverages)
    stats["frac_kmers_covered_gt1"] = (coverages > 1).sum().tolist() / len(coverages)
    stats["avg_kmer_coverage"] = coverages.mean().tolist()
    stats["median_kmer_coverage"] = int(np.median(coverages).tolist())
    stats["max_kmer_coverage"] = coverages.max().tolist()
    stats["frac_mult_weighted_kmers_covered"] = (
        multiplicities[coverages > 0].sum() / multiplicities.sum()
    ).tolist()
    stats["nterm_kmer_cov"] = ((coverages[nterms] > 0).sum() / nterms.sum()).tolist()
    stats["cterm_kmer_cov"] = ((coverages[cterms] > 0).sum() / cterms.sum()).tolist()
    stats["mult_hist"] = compute_int_hist(multiplicities)
    stats["cov_hist"] = compute_int_hist(coverages)

    return stats
