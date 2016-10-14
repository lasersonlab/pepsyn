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


from Bio.Seq import Seq

from pepsyn.error import PepsynError

def tile(seq, length, overlap):
    if length <= 0:
        raise ValueError('length must be a positive integer')
    if overlap < 0:
        raise ValueError('overlap must be a nonnegative integer')
    if overlap >= length:
        raise ValueError('length must be greater than overlap')
    if len(seq) < length:
        raise ValueError('length cannot be larger than the tiled seq')
    end = length
    while end < len(seq):
        start = end - length
        yield (start, end, seq[start:end])
        end += length - overlap
    if end == len(seq):
        start = end - length
        yield (start, end, seq[start:end])
    else:
        start = len(seq) - length
        end = len(seq)
        yield (start, end, seq[start:end])


def reverse_translate(seq, codon_sampler):
    """
    seq is a Bio.Seq.Seq
    codon_sampler is a CodonSampler
    """
    codons = []
    for aa in seq:
        codons.append(codon_sampler.sample_codon(aa))
    return sum(codons, Seq('', codon_sampler.nucleotide_alphabet))


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
