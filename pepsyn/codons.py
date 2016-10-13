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


import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Data.CodonTable import standard_dna_table


class CodonUsage(object):

    def __init__(self, freq):
        """freq is a dict of Seq('NNN', alphabet) -> float mappings"""
        # TODO: input verification:
        # all codons
        # sums to 1
        self.freq = freq
        self.nucleotide_alphabet = freq.keys().__iter__().__next__().alphabet


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
            c = Seq(codon, self.table.nucleotide_alphabet)
            self.aa2codons.setdefault(aa, []).append(c)

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
        usage is CodonUsage; None means uniform
        """
        super().__init__(table)
        if usage is None:
            raise ValueError('Must provide a CodonUsage object')
        if table.nucleotide_alphabet != usage.nucleotide_alphabet:
            raise ValueError('table and usage need to use the same nucleotide'
                             'Alphabet')
        # precalculate distribution (p) for each amino acid
        self.aa2p = {}
        if self.usage is not None:
            for aa in self.table.protein_alphabet.letters:
                unnormed = np.asarray(
                    [self.usage.freq[c] for c in self.aa2codons[aa]])
                self.aa2p[aa] = unnormed / unnormed.sum()

    def sample_codon(self, aa):
        """
        aa is str for single-letter IUPAC AA
        """
        i = np.random.choice(range(len(self.aa2codons[aa])), p=self.aa2p[aa])
        return self.aa2codons[aa][i]


# http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=37762
ecoli_codon_usage = CodonUsage({
    Seq('GGG', unambiguous_dna): 12.32,
    Seq('GGA', unambiguous_dna): 13.61,
    Seq('GGT', unambiguous_dna): 23.72,
    Seq('GGC', unambiguous_dna): 20.58,
    Seq('GAG', unambiguous_dna): 19.37,
    Seq('GAA', unambiguous_dna): 35.06,
    Seq('GAT', unambiguous_dna): 33.75,
    Seq('GAC', unambiguous_dna): 17.86,
    Seq('GTG', unambiguous_dna): 19.87,
    Seq('GTA', unambiguous_dna): 13.07,
    Seq('GTT', unambiguous_dna): 21.56,
    Seq('GTC', unambiguous_dna): 13.09,
    Seq('GCG', unambiguous_dna): 21.09,
    Seq('GCA', unambiguous_dna): 23.00,
    Seq('GCT', unambiguous_dna): 18.89,
    Seq('GCC', unambiguous_dna): 21.63,
    Seq('AGG', unambiguous_dna): 3.96,
    Seq('AGA', unambiguous_dna): 7.11,
    Seq('AGT', unambiguous_dna): 13.19,
    Seq('AGC', unambiguous_dna): 14.27,
    Seq('AAG', unambiguous_dna): 15.30,
    Seq('AAA', unambiguous_dna): 37.21,
    Seq('AAT', unambiguous_dna): 29.32,
    Seq('AAC', unambiguous_dna): 20.26,
    Seq('ATG', unambiguous_dna): 23.75,
    Seq('ATA', unambiguous_dna): 13.33,
    Seq('ATT', unambiguous_dna): 29.58,
    Seq('ATC', unambiguous_dna): 19.40,
    Seq('ACG', unambiguous_dna): 13.64,
    Seq('ACA', unambiguous_dna): 15.14,
    Seq('ACT', unambiguous_dna): 13.09,
    Seq('ACC', unambiguous_dna): 18.94,
    Seq('TGG', unambiguous_dna): 13.39,
    Seq('TGA', unambiguous_dna): 1.15,
    Seq('TGT', unambiguous_dna): 5.86,
    Seq('TGC', unambiguous_dna): 5.48,
    Seq('TAG', unambiguous_dna): 0.32,
    Seq('TAA', unambiguous_dna): 2.00,
    Seq('TAT', unambiguous_dna): 21.62,
    Seq('TAC', unambiguous_dna): 11.69,
    Seq('TTG', unambiguous_dna): 12.91,
    Seq('TTA', unambiguous_dna): 17.43,
    Seq('TTT', unambiguous_dna): 24.36,
    Seq('TTC', unambiguous_dna): 13.95,
    Seq('TCG', unambiguous_dna): 8.18,
    Seq('TCA', unambiguous_dna): 13.09,
    Seq('TCT', unambiguous_dna): 13.08,
    Seq('TCC', unambiguous_dna): 9.71,
    Seq('CGG', unambiguous_dna): 7.91,
    Seq('CGA', unambiguous_dna): 4.81,
    Seq('CGT', unambiguous_dna): 15.93,
    Seq('CGC', unambiguous_dna): 14.04,
    Seq('CAG', unambiguous_dna): 26.74,
    Seq('CAA', unambiguous_dna): 14.42,
    Seq('CAT', unambiguous_dna): 12.41,
    Seq('CAC', unambiguous_dna): 7.34,
    Seq('CTG', unambiguous_dna): 37.44,
    Seq('CTA', unambiguous_dna): 5.56,
    Seq('CTT', unambiguous_dna): 14.51,
    Seq('CTC', unambiguous_dna): 9.47,
    Seq('CCG', unambiguous_dna): 14.50,
    Seq('CCA', unambiguous_dna): 9.11,
    Seq('CCT', unambiguous_dna): 9.49,
    Seq('CCC', unambiguous_dna): 6.17})
