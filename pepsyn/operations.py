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


def assign_codons(seq, reverse_translator):
    """
    seq is a Bio.Seq.Seq
    reverse_translator is a ReverseTranslator
    """
    # TODO: check that protein alphabets match
    codons = []
    for aa in seq:
        codons.append(reverse_translator.sample_codon(aa))
    return sum(codons, Seq('', reverse_translator._nucleotide_alphabet))


def remove_cds_restriction_sites(seq, site, reverse_translator, cds_start=None,
                                 cds_end=None):
    """
    seq is a Bio.Seq.Seq
    site is a Bio.Seq.Seq
    reverse_translator is a ReverseTranslator
    start, end are int specifying the location of the CDS; defaults to full seq
    """
    site_start = seq.find(site)
    if site_start == -1:
        return seq
    if cds_start is None:
        cds_start = 0
    if cds_end is None:
        cds_end = len(seq)
    if (cds_end - cds_start) % 3 != 0:
        raise PepsynError('CDS length is not multiple of 3')
    if site_start + len(site) <= cds_start or site_start >= cds_end:
        raise PepsynError('Restriction site does not overlap CDS')
    mutseq = seq.tomutable()
    # computes offsets for the site boundaries to align them to coding frame
    start_offset = (site_start - cds_start) % 3
    end_offset = (site_start + len(site) - cds_start) % 3
    repl_start = max(site_start - start_offset, 0)
    repl_end = min(site_start + len(site) - end_offset + 3, len(seq))
    trans = seq[repl_start:repl_end].translate(table=reverse_translator._table)
    candidate = assign_codons(trans, reverse_translator)
    mutseq[repl_start:repl_end] = candidate
    if mutseq.find(site) == -1:
        return mutseq.toseq()
    else:
        # TODO: this can cause a maximum recursion error
        return remove_cds_restriction_sites(
            mutseq.toseq(), site, reverse_translator, cds_start=cds_start,
            cds_end=cds_end)
