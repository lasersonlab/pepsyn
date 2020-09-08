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

import warnings
from copy import copy
from itertools import product
from math import fsum

import numpy as np
from Bio.Data.CodonTable import standard_dna_table
from Bio.Seq import Seq

amber_codon = Seq("TAG")
ochre_codon = Seq("TAA")
opal_codon = Seq("TGA")


class CodonUsage(object):
    def __init__(self, weights):
        """weights is a dict of Seq('NNN', alphabet) -> float mappings

        weights does not need to be normalized, but must include all 64 codons
        """
        self.nucleotide_alphabet = standard_dna_table.nucleotide_alphabet

        # verify presence of all 64 codons
        all_64 = {Seq("".join(x)) for x in product(self.nucleotide_alphabet, repeat=3)}
        if all_64 != set(weights.keys()):
            raise ValueError("values must include all 64 codons")

        self.freq = copy(weights)
        self._renormalize_weights()

    def _renormalize_weights(self):
        total = fsum(self.freq.values())
        for codon in self.freq:
            self.freq[codon] = self.freq[codon] / total


# TODO: These fns seem like they need to be organized in a more principled
# manner
def zero_non_amber_stops(usage):
    """Returns CodonUsage that zeros-out non-amber stop codons"""
    # TODO: it doesn't take into account alphabet
    freq = usage.freq.copy()
    freq[ochre_codon] = 0.0
    freq[opal_codon] = 0.0
    return CodonUsage(freq)


def zero_low_freq_codons(usage, table, freq_threshold=0.01):
    """Returns CodonUsage that zeros low-freq codons unless the AA is elim"""
    freq = usage.freq.copy()
    common_codons = {str(c) for (c, f) in freq.items() if f >= freq_threshold}
    for aa in table.protein_alphabet + "*":
        curr_codons = {c for (c, a) in table.forward_table.items() if a == aa}
        if len(common_codons & curr_codons) == 0:
            continue
        rare_codons = curr_codons - common_codons
        for c in rare_codons:
            freq[Seq(c)] = 0.0
    return CodonUsage(freq)


class CodonSampler(object):
    def __init__(self, table=None):
        """
        table is Bio.Data.CodonTable.NCBICodonTableDNA
        """
        self.table = table if table is not None else standard_dna_table
        self.nucleotide_alphabet = self.table.nucleotide_alphabet
        self.protein_alphabet = self.table.protein_alphabet
        # generate aa -> codon mapping
        self.aa2codons = {}
        for (codon, aa) in self.table.forward_table.items():
            c = Seq(codon)
            self.aa2codons.setdefault(aa, []).append(c)
        for codon in self.table.stop_codons:
            c = Seq(codon)
            self.aa2codons.setdefault("*", []).append(c)

    def sample_codon(self, aa):
        """
        aa is str for single-letter IUPAC AA
        """
        raise NotImplementedError()


class UniformCodonSampler(CodonSampler):
    def __init__(self, table=None):
        """
        table is Bio.Data.CodonTable.NCBICodonTableDNA
        """
        super().__init__(table)

    def sample_codon(self, aa):
        i = np.random.randint(len(self.aa2codons[aa]))
        return self.aa2codons[aa][i]


class FreqWeightedCodonSampler(CodonSampler):
    def __init__(self, table=None, usage=None):
        """
        table is Bio.Data.CodonTable.NCBICodonTableDNA
        usage is CodonUsage
        """
        super().__init__(table)
        if usage is None:
            raise ValueError("Must provide a CodonUsage object")
        self.usage = usage
        if self.table.nucleotide_alphabet != self.usage.nucleotide_alphabet:
            raise ValueError("table and usage need to use the same nucleotide Alphabet")

        # precalculate amino acid distributions (incl stop codon)
        self.aa2p = {}
        for aa in self.table.protein_alphabet + "*":
            unnormed = np.asarray([self.usage.freq[c] for c in self.aa2codons[aa]])
            self.aa2p[aa] = unnormed / unnormed.sum()

    def sample_codon(self, aa):
        """
        aa is str for single-letter IUPAC AA
        """
        i = np.random.choice(range(len(self.aa2codons[aa])), p=self.aa2p[aa])
        return self.aa2codons[aa][i]


with warnings.catch_warnings():
    # to supress biopythons warnings on hashing Seq objects
    warnings.simplefilter("ignore")

    # http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=37762
    ecoli_codon_usage = CodonUsage(
        {
            Seq("GGG"): 12.32,
            Seq("GGA"): 13.61,
            Seq("GGT"): 23.72,
            Seq("GGC"): 20.58,
            Seq("GAG"): 19.37,
            Seq("GAA"): 35.06,
            Seq("GAT"): 33.75,
            Seq("GAC"): 17.86,
            Seq("GTG"): 19.87,
            Seq("GTA"): 13.07,
            Seq("GTT"): 21.56,
            Seq("GTC"): 13.09,
            Seq("GCG"): 21.09,
            Seq("GCA"): 23.00,
            Seq("GCT"): 18.89,
            Seq("GCC"): 21.63,
            Seq("AGG"): 3.96,
            Seq("AGA"): 7.11,
            Seq("AGT"): 13.19,
            Seq("AGC"): 14.27,
            Seq("AAG"): 15.30,
            Seq("AAA"): 37.21,
            Seq("AAT"): 29.32,
            Seq("AAC"): 20.26,
            Seq("ATG"): 23.75,
            Seq("ATA"): 13.33,
            Seq("ATT"): 29.58,
            Seq("ATC"): 19.40,
            Seq("ACG"): 13.64,
            Seq("ACA"): 15.14,
            Seq("ACT"): 13.09,
            Seq("ACC"): 18.94,
            Seq("TGG"): 13.39,
            Seq("TGA"): 1.15,
            Seq("TGT"): 5.86,
            Seq("TGC"): 5.48,
            Seq("TAG"): 0.32,
            Seq("TAA"): 2.00,
            Seq("TAT"): 21.62,
            Seq("TAC"): 11.69,
            Seq("TTG"): 12.91,
            Seq("TTA"): 17.43,
            Seq("TTT"): 24.36,
            Seq("TTC"): 13.95,
            Seq("TCG"): 8.18,
            Seq("TCA"): 13.09,
            Seq("TCT"): 13.08,
            Seq("TCC"): 9.71,
            Seq("CGG"): 7.91,
            Seq("CGA"): 4.81,
            Seq("CGT"): 15.93,
            Seq("CGC"): 14.04,
            Seq("CAG"): 26.74,
            Seq("CAA"): 14.42,
            Seq("CAT"): 12.41,
            Seq("CAC"): 7.34,
            Seq("CTG"): 37.44,
            Seq("CTA"): 5.56,
            Seq("CTT"): 14.51,
            Seq("CTC"): 9.47,
            Seq("CCG"): 14.50,
            Seq("CCA"): 9.11,
            Seq("CCT"): 9.49,
            Seq("CCC"): 6.17,
        }
    )
