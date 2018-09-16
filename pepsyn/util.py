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

import collections
import itertools

import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import _verify_alphabet
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Restriction.Restriction import enzymedict


def site2dna(site):
    """Convert "site" to a DNA sequence

    site is a str
    returns a Bio.Seq.Seq

    First tries to match site to a restriction enzyme. On failure, tries to
    convert to DNA sequence and checks strict alphabet
    """
    if site in enzymedict:
        dna = Seq(enzymedict[site]["site"], unambiguous_dna)
    else:
        dna = Seq(site, unambiguous_dna)
    if not _verify_alphabet(dna):
        raise ValueError("site is not recognized enzyme and not strict DNA")
    return dna


def sliding_window(n, seq):
    """ A sequence of overlapping subsequences

    Stolen pytoolz impl
    """
    return zip(
        *(
            collections.deque(itertools.islice(it, i), 0) or it
            for i, it in enumerate(itertools.tee(seq, n))
        )
    )


# fmt: off
# https://github.com/lh3/readfq
def readfq(fp): # this is a generator function  # pragma: no cover
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break
# fmt: on


def compute_int_hist(values):
    (hist, bin_edges) = np.histogram(values, bins=range(max(values) + 1))
    return {"hist": hist.tolist(), "bin_edges": bin_edges.tolist()}


def compute_float_hist(values):
    (hist, bin_edges) = np.histogram(values, bins="auto")
    return {"hist": hist.tolist(), "bin_edges": bin_edges.tolist()}
