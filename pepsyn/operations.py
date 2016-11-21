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


from Bio.Seq import Seq
from Bio.Data.IUPACData import protein_letters


from pepsyn.error import PepsynError


ambiguous_protein_values = {
    'B': 'DN',
    'X': protein_letters,
    'Z': 'EQ',
    'J': 'LI',
    'U': 'C',  # selenocysteine
    'O': 'K',  # pyrrolysine
}
extended_protein_letters = ''.join(ambiguous_protein_values.keys())


def tile(seq, length, overlap):
    """Generator of tiles with specified length across a sequence

    Generates the longest possible sequence of tiles, but may leave non-covered
    segments at the end of the sequence.  Sequence must be at least as long as
    tile length

    seq is a Bio.Seq.Seq
    length, overlap are int
    """
    if length <= 0:
        raise ValueError('length must be a positive integer')
    if overlap < 0:
        raise ValueError('overlap must be a nonnegative integer')
    if overlap >= length:
        raise ValueError('length must be greater than overlap')
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
        seq += '*'
    return seq[-length:]


def reverse_translate(seq, codon_sampler):
    """
    seq is a Bio.Seq.Seq
    codon_sampler is a CodonSampler
    """
    codons = []
    for aa in seq:
        codons.append(codon_sampler.sample_codon(aa))
    return sum(codons, Seq('', codon_sampler.nucleotide_alphabet))


def _ggsg_generator():
    """generates infinite sequence of 'G', 'G', 'S', 'G'"""
    for aa in cycle('GGSG'):
        yield aa


def x_to_ggsg(seq):
    """replace Xs with a Serine-Glycine linker (GGSG pattern)

    seq and return value are Bio.Seq.Seq
    """
    if 'X' not in seq:
        return seq
    replacement = []
    ggsg = _ggsg_generator()
    for aa in seq:
        if aa != 'X':
            replacement.append(aa)
            # restart linker iterator for next stretch of Xs
            ggsg = _ggsg_generator()
        else:
            replacement.append(next(ggsg))
    return Seq(''.join(replacement), seq.alphabet)


def pad_ggsg(seq, length, terminus='C'):
    """pad seq with Serine-Glycine linker (GGSG pattern)

    seq and return are Bio.Seq.Seq
    """
    if len(seq) >= length:
        return seq
    ggsg = _ggsg_generator()
    pad = ''.join([next(ggsg) for _ in range(length - len(seq))])
    if terminus == 'C':
        return seq + pad
    elif terminus == 'N':
        return pad + seq
    else:
        raise ValueError('terminus must be "N" or "C"')


def disambiguate_iupac_aa(seq):
    """generator
    seq is Bio.Seq.Seq
    """
    match = re.search('[{}]'.format(extended_protein_letters), seq.tostring())
    if match is None:
        yield seq
    else:
        idx = match.start()
        for aa in ambiguous_protein_values[seq[idx]]:
            disambiguated = seq[:idx] + aa + seq[idx + 1:]
            for s in disambiguate_iupac_aa(disambiguated):
                yield s


def remove_site_from_cds(seq, site, codon_sampler, cds_start=None,
                         cds_end=None):
    """
    seq is a Bio.Seq.Seq
    site is a Bio.Seq.Seq
    codon_sampler is a CodonSampler
    cds_start, cds_end are int; defaults to full seq
    """
    site_start = seq.find(site)
    if site_start == -1:
        return seq
    site_end = site_start + len(site)
    if cds_start is None:
        cds_start = 0
    if cds_end is None:
        cds_end = len(seq)
    if site_end <= cds_start or site_start >= cds_end:
        # site doesn't overlap the CDS => no-op
        return seq
    if (cds_end - cds_start) % 3 != 0:
        raise PepsynError('CDS length is not multiple of 3')

    # computes offsets for the site boundaries to align them to coding frame
    start_offset = (site_start - cds_start) % 3
    # negative because the end coord must be pushed to the right
    end_offset = -(site_end - cds_start) % 3
    repl_start = max(site_start - start_offset, cds_start)
    repl_end = min(site_end + end_offset, cds_end)

    # compute candidate recoded sequence
    trans = seq[repl_start:repl_end].translate(table=codon_sampler.table)
    replacement = reverse_translate(trans, codon_sampler)
    candidate = seq[:repl_start] + replacement + seq[repl_end:]
    if candidate.find(site) == -1:
        return candidate
    else:
        # TODO: this can cause a maximum recursion error
        return remove_site_from_cds(
            candidate, site, codon_sampler, cds_start=cds_start,
            cds_end=cds_end)
