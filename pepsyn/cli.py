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


from click import command, option, Path

from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Restriction.Restriction import enzymedict

@command(context_settings={'help_option_names': ['-h', '--help']})
@option('--input', '-i', type=Path(exists=True, allow_dash=True),
        help='Input protein sequences (fasta)')
@option('--output', '-o', type=Path(allow_dash=True),
        help='Output DNA oligos (fasta)')
@option('--length', '-l', type=int, help='Length of peptides (amino acids)')
@option('--overlap', '-p', type=int, help='Overlap of peptides (amino acids)')
@option('--prefix', default='', help='DNA sequence to prefix each oligo')
@option('--suffix', default='', help='DNA sequence to suffix each oligo')
@option('--remove', '-r', multiple=True, help=(
    'Remove restriction site. Use common name (e.g., EcoRI) or DNA sequence. '
    'Can be specified multiple times for multiple enzymes.'))
#TODO: take options for codon usage table/genetic code
def cli(input, output, length, overlap, prefix, suffix, remove):
    """pepsyn - peptide synthesis design"""
    if input is None:
        raise ValueError('Please specify an input file')
    if output is None:
        raise ValueError('Please specify an output file')
    prefix = Seq(prefix, unambiguous_dna)
    suffix = Seq(prefix, unambiguous_dna)
    remove_sites = []
    for site in raw_sites:
        if site in enzymedict:
            remove_sites.append(Seq(enzymedict[site]['site'], unambiguous_dna))
        else:
            remove_sites.append(Seq(site, unambiguous_dna))

    reverse_translator = ReverseTranslator(ecoli_codon_usage)

    with open(output_path, 'w') as op:
        for sr in SeqIO.parse(input_path, 'fasta'):
            for (s, e, t) in tile(sr.seq, length, overlap):
                dna = assign_codons(t, reverse_translator)
                dna = prefix + dna + suffix
                for site in remove_sites:
                    dna = remove_cds_restriction_sites(
                        dna, site, reverse_translator, len(prefix),
                        len(prefix) + len(dna))
                output_record = SeqRecord(dna, '{}|{}-{}'.format(id, s, e))
                SeqIO.write([output_record], op, 'fasta')
