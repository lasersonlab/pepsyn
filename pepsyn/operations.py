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

"""Operations on sequences.

This module should generally operate on Bio.Seq.Seq objects.
"""

import re
from itertools import cycle

from Bio.Data.IUPACData import (
    ambiguous_dna_values,
    protein_letters,
    unambiguous_dna_letters,
)
from Bio.Seq import Seq
from pygtrie import CharTrie

from pepsyn.error import PepsynError
from pepsyn.util import compute_float_hist, compute_int_hist

ambiguous_protein_values = {
    "B": "DN",
    "X": protein_letters,
    "Z": "EQ",
    "J": "LI",
    "U": "C",  # selenocysteine
    "O": "K",  # pyrrolysine
}
extended_protein_letters = "".join(ambiguous_protein_values.keys())
extended_dna_letters = set(ambiguous_dna_values.keys()) - set(unambiguous_dna_letters)


def tile(seq, length, overlap):
    """Generator of tiles with specified length across a sequence

    Generates the longest possible sequence of tiles, but may leave non-covered
    segments at the end of the sequence.  Sequence must be at least as long as
    tile length

    seq is a Bio.Seq.Seq
    length, overlap are int
    """
    if length <= 0:
        raise ValueError("length must be a positive integer")
    if overlap < 0:
        raise ValueError("overlap must be a nonnegative integer")
    if overlap >= length:
        raise ValueError("length must be greater than overlap")
    seqlen = len(seq)
    end = length
    while end < seqlen:
        start = end - length
        yield (start, end, seq[start:end])
        end += length - overlap
    if end == seqlen:
        start = end - length
        yield (start, end, seq[start:end])


def ctermpep(seq, length, add_stop=False):
    """Get C-terminal peptide

    add_stop will add a '*' to the peptide (total length still equal to
    `length`)

    If length is bigger than seq, it will return the entire seq

    seq is a Bio.Seq.Seq
    """
    if add_stop:
        seq += "*"
    return seq[-length:]


def reverse_translate(seq, codon_sampler):
    """
    seq is a Bio.Seq.Seq
    codon_sampler is a CodonSampler
    """
    codons = []
    for aa in seq:
        codons.append(codon_sampler.sample_codon(aa))
    return sum(codons, Seq(""))


def _ggsg_generator():
    """generates infinite sequence of 'G', 'G', 'S', 'G'"""
    for aa in cycle("GGSG"):
        yield aa


def x_to_ggsg(seq):
    """replace Xs with a Serine-Glycine linker (GGSG pattern)

    seq and return value are strings
    """
    if "X" not in seq:
        return seq
    replacement = []
    ggsg = _ggsg_generator()
    for aa in seq:
        if aa != "X":
            replacement.append(aa)
            # restart linker iterator for next stretch of Xs
            ggsg = _ggsg_generator()
        else:
            replacement.append(next(ggsg))
    return "".join(replacement)


def pad_ggsg(seq, length, terminus="C"):
    """pad seq with Serine-Glycine linker (GGSG pattern)

    seq and return are Bio.Seq.Seq
    """
    if len(seq) >= length:
        return seq
    ggsg = _ggsg_generator()
    pad = "".join([next(ggsg) for _ in range(length - len(seq))])
    if terminus == "C":
        return seq + pad
    elif terminus == "N":
        return pad + seq
    else:
        raise ValueError('terminus must be "N" or "C"')


def disambiguate_iupac_string(seq, ambiguous_letters, disambiguation):
    """generator
    seq is Bio.Seq.Seq

    ambiguous_letters is string containing ambiguous IUPAC codes

    disambiguation is dict with key from ambiguous_letters and values the
    strings containing the unambiguous versions (e.g.,
    Bio.Data.IUPACData.ambiguous_dna_values)

    """
    match = re.search("[{}]".format(ambiguous_letters), str(seq))
    if match is None:
        yield seq
    else:
        idx = match.start()
        for letter in disambiguation[seq[idx]]:
            disambiguated = seq[:idx] + letter + seq[idx + 1 :]
            for s in disambiguate_iupac_string(
                disambiguated, ambiguous_letters, disambiguation
            ):
                yield s


def disambiguate_iupac_dna(seq):
    """generator
    seq is Bio.Seq.Seq
    """
    for s in disambiguate_iupac_string(seq, extended_dna_letters, ambiguous_dna_values):
        yield s


def disambiguate_iupac_aa(seq):
    """generator
    seq is Bio.Seq.Seq
    """
    for s in disambiguate_iupac_string(
        seq, extended_protein_letters, ambiguous_protein_values
    ):
        yield s


def num_disambiguated_iupac_strings(seq, ambiguous_letters, disambiguation):
    """
    seq is Bio.Seq.Seq or str

    ambiguous_letters is string containing ambiguous IUPAC codes

    disambiguation is dict with key from ambiguous_letters and values the
    strings containing the unambiguous versions (e.g.,
    Bio.Data.IUPACData.ambiguous_dna_values)
    """
    n = 1
    for letter in str(seq):
        if letter in ambiguous_letters:
            n *= len(disambiguation[letter])
    return n


def num_disambiguated_iupac_dna(seq):
    return num_disambiguated_iupac_strings(
        seq, extended_dna_letters, ambiguous_dna_values
    )


def num_disambiguated_iupac_aa(seq):
    return num_disambiguated_iupac_strings(
        seq, extended_protein_letters, ambiguous_protein_values
    )


def recode(seq, codon_sampler):
    """
    seq is a Bio.Seq.Seq
    codon_sampler is a CodonSampler
    """
    translation = seq.translate(table=codon_sampler.table)
    recoded = reverse_translate(translation, codon_sampler)
    return recoded


def recode_sites_from_cds(
    seq, sites, codon_sampler, cds_start=None, cds_end=None, search_start=None
):
    """
    seq is a Bio.Seq.Seq
    site is a list[Bio.Seq.Seq]
    codon_sampler is a CodonSampler
    cds_start, cds_end are int; defaults to full seq; defines recodable region
    search_start is int; default full seq; where to start search for sites
    """
    # validate input
    search_start = search_start if search_start is not None else 0
    cds_start = cds_start if cds_start is not None else 0
    cds_end = cds_end if cds_end is not None else len(seq)
    if (cds_end - cds_start) % 3 != 0:
        raise PepsynError("CDS length is not multiple of 3")

    # initial search for site; handle recursion breaking cases here
    # search for all potential sites
    search_results = [(seq.find(site, search_start), len(site)) for site in sites]
    # filter out sites that don't match
    positive_matches = [tup for tup in search_results if tup[0] >= 0]
    if len(positive_matches) == 0:
        return seq
    # choose earliest matching site coord
    site_start = min([tup[0] for tup in positive_matches])
    # for any site that starts at site_start (could be multiple), choose the
    # largest site length
    site_end = site_start + max(
        [tup[1] for tup in positive_matches if tup[0] == site_start]
    )
    if site_end <= cds_start or site_start >= cds_end:
        # sites don't overlap the CDS => move on further down the seq
        return recode_sites_from_cds(
            seq,
            sites,
            codon_sampler,
            cds_start=cds_start,
            cds_end=cds_end,
            search_start=(site_start + 1),
        )
    # site needs to be recoded
    # computes offsets for the site boundaries to align them to coding frame
    start_offset = (site_start - cds_start) % 3
    # negative because the end coord must be pushed to the right
    end_offset = -(site_end - cds_start) % 3
    recode_start = max(site_start - start_offset, cds_start)
    recode_end = min(site_end + end_offset, cds_end)
    # compute candidate recoded sequence
    recoded_chunk = recode(seq[recode_start:recode_end], codon_sampler)
    candidate = seq[:recode_start] + recoded_chunk + seq[recode_end:]
    # TODO: this can cause a maximum recursion error if it never finds a
    # recoded sequence without the site recursively
    return recode_sites_from_cds(
        candidate,
        sites,
        codon_sampler,
        cds_start=cds_start,
        cds_end=cds_end,
        search_start=search_start,
    )


def recode_site_from_cds(
    seq, site, codon_sampler, cds_start=None, cds_end=None, search_start=None
):
    return recode_sites_from_cds(
        seq,
        [site],
        codon_sampler,
        cds_start=cds_start,
        cds_end=cds_end,
        search_start=search_start,
    )


def orf_stats(orfs):
    """compute orf stats

    orfs is name->seq dicts
    """
    import numpy as np

    orf_lens = np.asarray([len(o) for o in orfs.values()])
    ambiguity_factors = {n: num_disambiguated_iupac_aa(s) for (n, s) in orfs.items()}
    stats = {}
    stats["num_orfs"] = len(orfs)
    stats["total_orf_residues"] = orf_lens.sum().tolist()
    stats["orf_lens_hist"] = compute_int_hist(orf_lens)
    stats["avg_orf_len"] = orf_lens.mean().tolist()
    stats["max_orf_len"] = orf_lens.max().tolist()
    stats["min_orf_len"] = orf_lens.min().tolist()
    stats["num_orfs_internal_stops"] = sum(
        ["*" in s.rstrip("*") for s in orfs.values()]
    )
    stats["num_orfs_Xs"] = sum(["X" in s.upper() for s in orfs.values()])
    stats["max_ambig_factor"] = max(ambiguity_factors.values())
    stats["ambig_factor_hist"] = compute_int_hist(list(ambiguity_factors.values()))
    stats["top_5_ambig"] = list(
        map(
            list,
            sorted(ambiguity_factors.items(), key=lambda tup: tup[1], reverse=True)[:5],
        )
    )

    return stats


def tile_stats(orfs, tiles):
    """compute tile stats

    orfs and tiles are name->seq dicts

    NOTE: for prefix trie stats (e.g., num of tiles per orf), it is assumed the
    orf name is a prefix to the name of a tile from that orf
    """
    import numpy as np

    tile_lens = np.asarray([len(t) for t in tiles.values()])
    orf_lens = np.asarray([len(o) for o in orfs.values()])
    tile_size = int(round(np.median(tile_lens)).tolist())

    # compute tile counts for each orf
    orf_prefixes = CharTrie()
    for name in orfs:
        orf_prefixes[name] = True
    # ensure that no ORF name is a prefix for a different valid ORF
    for name in orfs:
        if len(orf_prefixes.keys(name)) != 1:
            print(orf_prefixes.keys(name))
            raise ValueError("some ORF name is a prefix for a different valid ORF")
    tile_prefixes = CharTrie()
    for name in tiles:
        tile_prefixes[name] = True
    # compute orf coverages
    orf_coverages = {}
    for (orf, seq) in orfs.items():
        orf_residues = len(seq)
        tile_residues = 0.0
        if tile_prefixes.has_subtrie(orf) or (orf in tile_prefixes):
            for tile in tile_prefixes.keys(orf):
                tile_residues += len(tiles[tile])
        orf_coverages[orf] = tile_residues / orf_residues

    stats = {}
    stats["tile_size"] = tile_size
    stats["num_tiles"] = len(tiles)
    stats["total_tile_residues"] = tile_lens.sum().tolist()
    stats["avg_orf_coverage"] = tile_lens.sum().tolist() / orf_lens.sum().tolist()
    stats["num_orfs_smaller_than_tile_size"] = (orf_lens < tile_size).sum().tolist()
    stats["approx_num_tiles_naive_1x_tiling"] = (
        np.ceil(orf_lens / tile_size).sum().tolist()
    )
    stats["avg_orf_coverage"] = sum(orf_coverages.values()) / len(orf_coverages)
    stats["max_tiles_per_len_normed_orf"] = max(orf_coverages.values())
    stats["tile_len_hist"] = compute_int_hist(tile_lens)
    # what is the tile coverage of each ORF (tot tile residues / orf residues)
    # tiles are assigned to ORFs if they share a name
    stats["orf_coverage_hist"] = compute_float_hist(list(orf_coverages.values()))
    stats["top_5_orf_cov"] = list(
        map(
            list,
            sorted(orf_coverages.items(), key=lambda tup: tup[1], reverse=True)[:5],
        )
    )
    stats["bot_5_orf_cov"] = list(
        map(list, sorted(orf_coverages.items(), key=lambda tup: tup[1])[:5])
    )

    return stats
