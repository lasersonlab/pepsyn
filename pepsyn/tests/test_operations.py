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

from pytest import raises
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import protein
import numpy as np

from pepsyn.operations import tile, assign_codons
from pepsyn.codons import ReverseTranslator, ecoli_codon_usage


protein_seq = Seq('METMSDYSKEVSEALSALRGELSALSAAISNTVRAGSYSAPVAKDCKAGHCDSKAVL',
                  protein)
short_protein_seq = Seq('METMSD', protein)
reverse_translator = ReverseTranslator(ecoli_codon_usage)


class TestTile(object):

    def test_nonoverlapping_perfect(self):
        length = 19
        overlap = 0
        tiles = [t[2] for t in tile(protein_seq, length, overlap)]
        assert len(tiles) == 3
        assert all([len(t) == length for t in tiles])
        assert tiles[0] == Seq('METMSDYSKEVSEALSALR', protein)
        assert tiles[1] == Seq('GELSALSAAISNTVRAGSY', protein)
        assert tiles[2] == Seq('SAPVAKDCKAGHCDSKAVL', protein)
        assert sum(tiles, Seq('', protein)) == protein_seq

    def test_nonoverlapping_imperfect(self):
        length = 25
        overlap = 0
        tiles = [t[2] for t in tile(protein_seq, length, overlap)]
        assert len(tiles) == 3
        assert all([len(t) == length for t in tiles])
        assert tiles[0] == Seq('METMSDYSKEVSEALSALRGELSAL', protein)
        assert tiles[1] == Seq('SAAISNTVRAGSYSAPVAKDCKAGH', protein)
        assert tiles[2] == Seq('VRAGSYSAPVAKDCKAGHCDSKAVL', protein)

    def test_overlapping_perfect(self):
        length = 21
        overlap = 3
        tiles = [t[2] for t in tile(protein_seq, length, overlap)]
        assert len(tiles) == 3
        assert all([len(t) == length for t in tiles])
        assert tiles[0] == Seq('METMSDYSKEVSEALSALRGE', protein)
        assert tiles[1] == Seq('RGELSALSAAISNTVRAGSYS', protein)
        assert tiles[2] == Seq('SYSAPVAKDCKAGHCDSKAVL', protein)

    def test_overlapping_imperfect(self):
        length = 22
        overlap = 3
        tiles = [t[2] for t in tile(protein_seq, length, overlap)]
        assert len(tiles) == 3
        assert all([len(t) == length for t in tiles])
        assert tiles[0] == Seq('METMSDYSKEVSEALSALRGEL', protein)
        assert tiles[1] == Seq('GELSALSAAISNTVRAGSYSAP', protein)
        assert tiles[2] == Seq('GSYSAPVAKDCKAGHCDSKAVL', protein)

    def test_length_longer_than_seq(self):
        length = len(protein_seq) + 5
        overlap = 5
        with raises(ValueError):
            tiles = list(tile(protein_seq, length, overlap))

    def test_overlap_longer_than_length(self):
        length = 10
        overlap = 15
        with raises(ValueError):
            tiles = list(tile(protein_seq, length, overlap))

    def test_negative_length(self):
        length = -10
        overlap = 5
        with raises(ValueError):
            tiles = list(tile(protein_seq, length, overlap))

    def test_zero_length(self):
        length = 0
        overlap = 5
        with raises(ValueError):
            tiles = list(tile(protein_seq, length, overlap))

    def test_negative_overlap(self):
        length = 20
        overlap = -3
        with raises(ValueError):
            tiles = list(tile(protein_seq, length, overlap))
