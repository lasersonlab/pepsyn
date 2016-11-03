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
        dna = Seq(enzymedict[site]['site'], unambiguous_dna)
    else:
        dna = Seq(site, unambiguous_dna)
    if not _verify_alphabet(dna):
        raise ValueError('site is not recognized enzyme and not strict DNA')
    return dna
