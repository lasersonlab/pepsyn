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


from click import group, command, option, argument, File, Choice

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Restriction.Restriction import enzymedict

from pepsyn.operations import (
    reverse_translate, remove_site_from_cds, tile as tile_op)
from pepsyn.codons import (
    FreqWeightedCodonSampler, UniformCodonSampler, ecoli_codon_usage)


@group(context_settings={'help_option_names': ['-h', '--help']})
def cli():
    """pepsyn -- peptide synthesis design"""
    pass


# reusable args
argument_input = argument('input', type=File('r'))
argument_output = argument('output', type=File('w'))


@cli.command()
@argument_input
@argument_output
@option('--length', '-l', type=int, help='Length of output oligos')
@option('--overlap', '-p', type=int, help='Overlap of oligos')
def tile(input, output, length, overlap):
    """tile a set of sequences"""
    for seqrecord in SeqIO.parse(input, 'fasta'):
        for (start, end, t) in tile_op(seqrecord.seq, length, overlap):
                output_title = '{}|{}-{}'.format(seqrecord.id, start, end)
                output_record = SeqRecord(t, output_title, description='')
                SeqIO.write(output_record, output, 'fasta')


@cli.command()
@argument_input
@argument_output
@option('--prefix', '-p', help='DNA sequence to prefix each oligo')
def prefix(input, output, prefix):
    """add a prefix to each sequence"""
    for seqrecord in SeqIO.parse(input, 'fasta'):
        newseq = prefix + seqrecord.seq
        seqrecord.seq = newseq
        SeqIO.write(seqrecord, output, 'fasta')


@cli.command()
@argument_input
@argument_output
@option('--suffix', '-s', help='DNA sequence to suffix each oligo')
def suffix(input, output, suffix):
    """add a suffix to each sequence"""
    for seqrecord in SeqIO.parse(input, 'fasta'):
        newseq = seqrecord.seq + suffix
        seqrecord.seq = newseq
        SeqIO.write(seqrecord, output, 'fasta')


@cli.command()
@argument_input
@argument_output
@option('--codon-table', '-t', default='standard',
        help='ONLY STANDARD TABLE IMPLEMENTED')
@option('--codon-usage', '-u', default='ecoli', help='ONLY ECOLI IMPLEMENTED')
@option('--sampler', default='weighted', show_default=True,
        type=Choice(['weighted', 'uniform']), help='Codon sampling method')
def revtrans(input, output, codon_table, codon_usage, sampler):
    """reverse translate amino acid sequences into DNA"""
    if sampler == 'weighted':
        codon_sampler = FreqWeightedCodonSampler(usage=ecoli_codon_usage)
    elif sampler == 'uniform':
        codon_sampler = UniformCodonSampler()
    for seqrecord in SeqIO.parse(input, 'fasta'):
        dna_id = seqrecord.id
        dna_seq = reverse_translate(seqrecord.seq, codon_sampler)
        SeqIO.write(SeqRecord(dna_seq, dna_id, description=''), output, 'fasta')


@cli.command()
@argument_input
@argument_output
@option('--site', help='Site to remove (e.g., EcoRI, AGCCT); case sensitive')
@option('--start', type=int, help='Start position of CDS (Pythonic coords)')
@option('--end', type=int, help='End position of CDS (Pythonic coords)')
@option('--codon-table', '-t', default='standard',
        help='ONLY STANDARD TABLE IMPLEMENTED')
@option('--codon-usage', '-u', default='ecoli', help='ONLY ECOLI IMPLEMENTED')
@option('--sampler', default='weighted', show_default=True,
        type=Choice(['weighted', 'uniform']), help='Codon sampling method')
def removesite(input, output, site, start, end, codon_table, codon_usage,
               sampler):
    """remove site from each sequence's CDS by recoding"""
    if sampler == 'weighted':
        codon_sampler = FreqWeightedCodonSampler(usage=ecoli_codon_usage)
    elif sampler == 'uniform':
        codon_sampler = UniformCodonSampler()

    if site in enzymedict:
        site = Seq(enzymedict[site]['site'], unambiguous_dna)
    else:
        site = Seq(site, unambiguous_dna)

    for seqrecord in SeqIO.parse(input, 'fasta'):
        id_ = seqrecord.id
        seq = remove_site_from_cds(seqrecord.seq, site, codon_sampler, start,
                                   end)
        SeqIO.write(SeqRecord(seq, id_, description=''), output, 'fasta')


@cli.command()
@argument_input
def stats(input):
    """compute some sequence statistics"""
    pass
