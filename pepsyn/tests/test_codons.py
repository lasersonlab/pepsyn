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

from math import isclose

from Bio.Data.CodonTable import standard_dna_table

from pepsyn.codons import (
    zero_non_amber_stops, zero_low_freq_codons, ecoli_codon_usage, amber_codon,
    ochre_codon, opal_codon)


class TestUsageManipulation(object):

    def test_zero_non_amber(self):
        zeroed_weight = (
            ecoli_codon_usage.freq[ochre_codon] +
            ecoli_codon_usage.freq[opal_codon])
        new_usage = zero_non_amber_stops(ecoli_codon_usage)
        for codon in new_usage.freq:
            if codon == ochre_codon or codon == opal_codon:
                continue
            inflated_freq = ecoli_codon_usage.freq[codon] / (1 - zeroed_weight)
            new_freq = new_usage.freq[codon]
            assert isclose(new_freq, inflated_freq)
        assert new_usage.freq[ochre_codon] == 0
        assert new_usage.freq[opal_codon] == 0
