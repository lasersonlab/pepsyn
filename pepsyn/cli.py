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

import os
import sys
from collections import Counter
from logging import captureWarnings
from math import ceil, floor, inf, log10
from os.path import join as pjoin

import yaml
from Bio import SeqIO
from Bio.Data.CodonTable import standard_dna_table
from Bio.SeqRecord import SeqRecord
from click import (
    Abort,
    Choice,
    File,
    Path,
    UsageError,
    argument,
    command,
    group,
    option,
    version_option,
)
from tqdm import tqdm, trange

from pepsyn import __version__
from pepsyn.codons import (
    FreqWeightedCodonSampler,
    UniformCodonSampler,
    ecoli_codon_usage,
    zero_low_freq_codons,
    zero_non_amber_stops,
)
from pepsyn.operations import ctermpep as cterm_oligo
from pepsyn.operations import (
    disambiguate_iupac_aa,
    num_disambiguated_iupac_aa,
    orf_stats,
    pad_ggsg,
    recode_site_from_cds,
    recode_sites_from_cds,
    reverse_translate,
)
from pepsyn.operations import tile as tile_op
from pepsyn.operations import tile_stats, x_to_ggsg
from pepsyn.util import readfq, site2dna, sliding_window

# biopython has a bunch of annoying warnings bc Seq comparisons changed
captureWarnings(True)


def print_fasta(sr, out):
    print(">{}\n{}".format(sr.id, str(sr.seq)), file=out)


@group(context_settings={"help_option_names": ["-h", "--help"]})
@version_option(__version__)
def cli():
    """peptide synthesis design"""
    pass


# reusable args
argument_input = argument("input", type=File("r"))
argument_output = argument("output", type=File("w"))


@cli.command(short_help="tile across input sequences")
@argument_input
@argument_output
@option("--length", "-l", type=int, help="Length of output oligos")
@option("--overlap", "-p", type=int, help="Overlap of oligos")
def tile(input, output, length, overlap):
    """Generate short (overlapping) sequence tiles from input sequences.

    Each sequence in the fasta input is converted into short tiles with given
    length and overlap and written out in fasta format.

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    Note: this tool drops "incomplete/short" last tiles if the length/overlap
    setting does not allow a tiling to perfectly cover a sequence. We recommend
    using ``ctermpep`` to explicitly capture the last tiles.

    """
    for (name, seq, qual) in tqdm(readfq(input), desc="tile", unit="seq"):
        for (start, end, t) in tile_op(seq, length, overlap):
            output_title = f"{name}|{start}-{end}"
            print(f">{output_title}\n{t}", file=output)


@cli.command(short_help="extract C-terminal peptide")
@argument_input
@argument_output
@option("--length", "-l", type=int, help="Length of output C-terminal oligo")
@option("--add-stop", "-s", is_flag=True, help="Add a stop codon to peptide")
def ctermpep(input, output, length, add_stop):
    """Extract the C-terminal peptide from each input sequence

    If an input sequence is shorter than the specified length, it will write
    out the entirety of the sequence.

    With the ``--add-stop`` option, an asterisk is appended to the input
    sequence and counts as one of the amino acids in terms of peptide length.
    For example, if requesting 56-aa peptides with a stop codon, the output
    will code for 55 amino acids and the stop.

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    for (name, seq, qual) in tqdm(readfq(input), desc="ctermpep", unit="seq"):
        oligo = cterm_oligo(seq, length, add_stop=add_stop)
        output_title = f"{name}|CTERM"
        if add_stop:
            output_title = f"{output_title}|STOP"
        print(f">{output_title}\n{oligo}", file=output)


@cli.command(short_help="strip stop codons from protein sequences")
@argument_input
@argument_output
def stripstop(input, output):
    """Strip stop "codons" from end of input protein sequences.

    Stop codons are assumed to be represented as "*".

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    for (name, seq, qual) in readfq(input):
        seq = seq.rstrip("*")
        print(f">{name}\n{seq}", file=output)


@cli.command(short_help="filter out protein sequences with stops")
@argument_input
@argument_output
def filterstop(input, output):
    """Filter out input sequences that contain stop codons (*).

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    for (name, seq, qual) in readfq(input):
        if "*" in seq:
            continue
        print(f">{name}\n{seq}", file=output)


@cli.command(short_help="size select sequences")
@argument_input
@argument_output
@option("-m", "--min-len", type=int, default=0, help="Min length of sequence to keep")
@option("-M", "--max-len", type=int, help="Max length of sequence to keep")
def filterlen(input, output, min_len, max_len):
    """Filter sequences of a given length.

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    if max_len is None:
        max_len = inf
    for (name, seq, qual) in readfq(input):
        if len(seq) < min_len:
            continue
        if len(seq) > max_len:
            continue
        print(f">{name}\n{seq}", file=output)


@cli.command(short_help="filter out duplicate sequences")
@argument_input
@argument_output
def uniq(input, output):
    """Filter out duplicate sequences.

    Only takes account of the sequences themselves. Arbitrarily picks one. This
    requires loading the entire file into RAM.

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    seqs = {}
    for (name, seq, qual) in readfq(input):
        if seq not in seqs:
            seqs[seq] = name
    for (seq, name) in seqs.items():
        print(f">{name}\n{seq}", file=output)


@cli.command(short_help="pad peptide to specified length")
@argument_input
@argument_output
@option("--length", "-l", type=int, help="Target length for peptide")
@option(
    "--n-term/--c-term", "-n/-c", "nterm", help="C- or N- terminus [default: C-term]"
)
def pad(input, output, length, nterm):
    """Pad protein sequence to a specified length by adding amino acids in the
    pattern of "GSGG".

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    terminus = "N" if nterm else "C"
    for (name, seq, qual) in readfq(input):
        padded = pad_ggsg(seq, length, terminus)
        pad_len = len(padded) - len(seq)
        if pad_len > 0:
            output_title = f"{name}|{terminus}-PADDED-{pad_len}"
        else:
            output_title = name
        print(f">{output_title}\n{padded}", file=output)


@cli.command(short_help="add a prefix to each sequence")
@argument_input
@argument_output
@option("--prefix", "-p", help="DNA sequence to prefix each oligo")
def prefix(input, output, prefix):
    """Add a prefix to each sequence

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    for (name, seq, qual) in readfq(input):
        newseq = prefix + seq
        print(f">{name}\n{newseq}", file=output)


@cli.command(short_help="add a suffix to each sequence")
@argument_input
@argument_output
@option("--suffix", "-s", help="DNA sequence to suffix each oligo")
def suffix(input, output, suffix):
    """Add a suffix to each sequence

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    for (name, seq, qual) in readfq(input):
        newseq = seq + suffix
        print(f">{name}\n{newseq}", file=output)


@cli.command(short_help="clip/truncate each sequence")
@argument_input
@argument_output
@option("--left", "-l", default=0, help="number of bases to clip from left")
@option("--right", "-r", default=0, help="number of bases to clip from right")
def clip(input, output, left, right):
    """Clip/truncate bases from the ends of each sequence

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    for (name, seq, qual) in readfq(input):
        stop = len(seq) - right
        print(f">{name}\n{seq[left:stop]}", file=output)


@cli.command(short_help="translate DNA to protein sequences")
@argument_input
@argument_output
@option("--codon-table", "-t", default="Standard", help="Specify the codon table")
@option(
    "--truncate-at-stop",
    "-x",
    is_flag=True,
    help="Truncate translation at first stop codon",
)
def translate(input, output, codon_table, truncate_at_stop):
    """Translate nucleotide sequences into protein

    Note: only the Standard codon table is currently implemented.

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.
    """
    for seqrecord in tqdm(SeqIO.parse(input, "fasta"), desc="translate", unit="seq"):
        aa_id = seqrecord.id
        aa_seq = seqrecord.seq.translate(table=codon_table)
        if truncate_at_stop:
            aa_seq = aa_seq.split("*")[0]
        print_fasta(SeqRecord(aa_seq, aa_id, description=""), output)


@cli.command(short_help="reverse translate protein to DNA sequences")
@argument_input
@argument_output
@option("--codon-table", "-t", default="Standard", help="Specify the codon table")
@option("--codon-usage", "-u", default="ecoli", help="Specify the codon usage")
@option(
    "--sampler",
    default="weighted",
    show_default=True,
    type=Choice(["weighted", "uniform"]),
    help="Codon sampling method",
)
@option(
    "--codon-freq-threshold", type=float, default=None, help="Minimum codon frequency"
)
@option("--amber-only", is_flag=True, help="Use only amber stop codon")
def revtrans(
    input, output, codon_table, codon_usage, sampler, codon_freq_threshold, amber_only
):
    """Reverse translate amino acid sequences into DNA

    This operation randomly samples codons for each amino acid, so multiple runs
    of this tool on the same input can produce different results

    Note: only the Standard codon table is currently implemented.
    Note: only the E. coli codon usage is currently implemented.

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    if sampler == "weighted":
        usage = ecoli_codon_usage
        if codon_freq_threshold is not None:
            # TODO: this is hardcoded in and there's a leaky abstraction here
            table = standard_dna_table
            usage = zero_low_freq_codons(usage, table, codon_freq_threshold)
        if amber_only:
            usage = zero_non_amber_stops(usage)
        codon_sampler = FreqWeightedCodonSampler(usage=usage)
    elif sampler == "uniform":
        codon_sampler = UniformCodonSampler()
    for seqrecord in tqdm(SeqIO.parse(input, "fasta"), desc="revtrans", unit="seq"):
        dna_id = seqrecord.id
        dna_seq = reverse_translate(seqrecord.seq, codon_sampler)
        print_fasta(SeqRecord(dna_seq, dna_id, description=""), output)


@cli.command(short_help="recode DNA sequence")
@argument_input
@argument_output
@option(
    "--site",
    multiple=True,
    help="Site to remove (e.g., EcoRI, AGCCT); case sens; allows multiple",
)
@option(
    "--clip-left",
    type=int,
    default=0,
    help="Number of bases to clip from start of sequence to get to CDS",
)
@option(
    "--clip-right",
    type=int,
    default=0,
    help="Number of bases to clip from end of sequence to get to CDS",
)
@option(
    "--codon-table", "-t", default="Standard", help="ONLY STANDARD TABLE IMPLEMENTED"
)
@option("--codon-usage", "-u", default="ecoli", help="ONLY ECOLI IMPLEMENTED")
@option(
    "--sampler",
    default="weighted",
    show_default=True,
    type=Choice(["weighted", "uniform"]),
    help="Codon sampling method",
)
@option(
    "--codon-freq-threshold", type=float, default=None, help="Minimum codon frequency"
)
@option("--amber-only", is_flag=True, help="Use only amber stop codon")
def recodesite(
    input,
    output,
    site,
    clip_left,
    clip_right,
    codon_table,
    codon_usage,
    sampler,
    codon_freq_threshold,
    amber_only,
):
    """Recode a DNA sequence to remove a particular site (e.g., restriction site)

    The site needs to be recognized by Biopython, or it will be treated as a DNA
    sequence. The clipping options should determine the boundaries of the coding
    sequence, which will correspond to the part of the sequence that is
    "recodable".

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    if sampler == "weighted":
        usage = ecoli_codon_usage
        if codon_freq_threshold is not None:
            # TODO: this is hardcoded in and there's a leaky abstraction here
            table = standard_dna_table
            usage = zero_low_freq_codons(usage, table, codon_freq_threshold)
        if amber_only:
            usage = zero_non_amber_stops(usage)
        codon_sampler = FreqWeightedCodonSampler(usage=usage)
    elif sampler == "uniform":
        codon_sampler = UniformCodonSampler()

    sites = [site2dna(s) for s in site]
    # sites is now a list[Bio.Seq.Seq]

    for seqrecord in SeqIO.parse(input, "fasta"):
        id_ = seqrecord.id
        cds_start = clip_left
        cds_end = len(seqrecord) - clip_right
        seq = recode_sites_from_cds(
            seqrecord.seq, sites, codon_sampler, cds_start, cds_end
        )
        print_fasta(SeqRecord(seq, id_, description=""), output)


@cli.command(short_help="replace Xs with linker")
@argument_input
@argument_output
def x2ggsg(input, output):
    """Replace stretches of Xs with Serine-Glycine linker (in a GGSG pattern)

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    for (name, seq, qual) in readfq(input):
        replacement = x_to_ggsg(seq)
        if replacement != seq:
            output_title = f"{name}|withGSlinker"
        else:
            output_title = name
        print(f">{output_title}\n{replacement}", file=output)


@cli.command(short_help="replace ambiguous AAs with unambiguous ones")
@argument_input
@argument_output
def disambiguateaa(input, output):
    """Replace IUPAC ambiguous amino acids with unambiguous ones

    Specifically, make the following replacements:
    B => DN
    X => ACDEFGHIKLMNPQRSTVWY
    Z => EQ
    J => LI,
    U => C (selenocysteine)
    O => K (pyrrolysine)

    If there are multiple possible replacements, this operation will output a
    sequence for each possible option. Use caution with sequences that are
    highly ambiguous (e.g., with many Xs), as in this case a single sequence
    could lead to an explosion in the output.

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    for (name, ambig, qual) in readfq(input):
        n = num_disambiguated_iupac_aa(ambig)
        digits = floor(log10(n)) + 1
        fmt = f"{name}|disambig_{{:0{digits}d}}"
        for (i, unambig) in enumerate(disambiguate_iupac_aa(ambig)):
            if n > 1:
                name = fmt.format(i + 1)
            print(f">{name}\n{unambig}", file=output)


@cli.command(short_help="find DNA site in sequences")
@argument_input
@option("--site", help="Site to find (e.g., EcoRI, AGCCT); case sensitive")
@option(
    "--clip-left",
    type=int,
    default=0,
    help="Number of bases to clip from start of sequence to get to CDS",
)
@option(
    "--clip-right",
    type=int,
    default=0,
    help="Number of bases to clip from end of sequence to get to CDS",
)
def findsite(input, site, clip_left, clip_right):
    """Find locations of a site in a DNA sequences

    If a sequence matches the specified site, write out its name and location.
    Used as a diagnostic to confirm that a particular DNA site (e.g.,
    restriction enzyme) is absent from a set of sequences. Because there may be
    adaptor sequences that contain such a site by design, the clipping option
    allows the search to be restricted. Note that a site is searched if it
    overlaps with the valid region even by one base (i.e., a site can match if
    it is mostly outside the clipped region, as long as it overlaps the target
    search region).

    INPUT is a path to fasta file or "-" to specify STDIN.

    """
    query = str(site2dna(site))
    for (name, seq, qual) in readfq(input):
        start = clip_left
        end = len(seq) - clip_right
        idx = seq[start:end].find(query)
        if idx >= 0:
            print(f"{name}|{site}|{idx + start}", file=sys.stdout)


@cli.command(short_help="build de Bruijn graph")
@argument_input
@argument("output_path", type=Path(exists=False))
@option("-k", "--kmer-size", type=int, required=True, help="k-mer size")
def builddbg(input, output_path, kmer_size):
    """Build a de bruijn graph on a set of protein sequences

    This process ignores input sequences shorter than the specified kmer size.

    If the output path ends with .gz, the output will be compressed.

    INPUT is a path to fasta file or "-" to specify STDIN. OUTPUT_PATH must
    point to a valid path.

    """
    try:
        import networkx as nx
        from pepsyn.dbg import fasta_handle_to_dbg
    except ImportError:
        raise Abort("builddbg requires NetworkX")

    with tqdm(desc="building dbg") as pbar:
        dbg = fasta_handle_to_dbg(input, kmer_size, tqdm=pbar, ignore_short=True)
    nx.write_gpickle(dbg, output_path)


@cli.command(short_help="sample tiles from de Bruijn graph")
@argument_input
@argument_output
@option("-t", "--tile-size", type=int, required=True, help="tile size")
@option(
    "-d",
    "--dbg-path",
    type=Path(exists=True, dir_okay=False),
    required=True,
    help="path to prebuilt de bruijn graph (pickle file; .gz okay)",
)
@option(
    "-c",
    "--kmer-cov",
    type=float,
    help="target avg k-mer coverage (incompatible with -n)",
)
@option(
    "-n", "--num-tiles", type=int, help="number of output tiles (incompatible with -c)"
)
@option(
    "-p",
    "--preselected-tiles-path",
    type=File("r"),
    help="fasta file containing previously selected tiles",
)
def greedykmercov(
    input, output, tile_size, dbg_path, kmer_cov, num_tiles, preselected_tiles_path
):
    """Select protein tiles from de Bruijn graph by maximizing k-mer coverage

    Each tile is a fragment of an observed input ORF. Either the total number of
    output tiles can be specified, or the average target k-mer coverage. If
    there is already a pre-selected set of tiles chosen through some other
    method, specifying them will initialize the de Bruijn graph to reflect the
    preexisting k-mer coverage.

    NOTE: ORFS shorter than tile-size are sampled, but ORFs shorter than
    kmer-size are ignored. (Use pepsyn filterlen to select short tiles.)

    Preselected tiles are given to this command to keep it from choosing the
    tile again. This does not affect the number of tiles output, however, nor
    their distribution across the DBG components when ignoring the preselected
    tiles.

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    # test input/context
    try:
        import networkx as nx
        from pepsyn.dbg import gen_kmers, setreduce_attr, sum_attr
    except ImportError:
        raise Abort("greedykmercov requires NetworkX")
    try:
        import numpy as np
    except ImportError:
        raise Abort("greedykmercov requires NumPy")
    if kmer_cov and num_tiles:
        raise UsageError("Set -c/--kmer-cov OR -n/--num-tiles but not both")
    if not kmer_cov and not num_tiles:
        raise UsageError("Must set one of -c/--kmer-cov OR -n/--num-tiles")

    # load orfs
    orfs = {name: seq for (name, seq, qual) in readfq(input)}

    # load dbg
    dbg = nx.read_gpickle(dbg_path)
    kmer_size = len(next(iter(dbg)))
    if kmer_size > tile_size:
        raise UsageError("kmer-size > tile_size")
    kmers_remaining = len(dbg)
    num_components = nx.number_weakly_connected_components(dbg)
    if num_tiles:
        tiles_remaining = num_tiles

    # load preselected tiles
    if preselected_tiles_path:
        preselected_tiles = [
            seq for (name, seq, qual) in readfq(preselected_tiles_path)
        ]
        preselected_kmer_counts = Counter(
            [
                kmer
                for tile in preselected_tiles
                for kmer in gen_kmers(tile, kmer_size, yield_short=True)
            ]
        )
    else:
        preselected_tiles = []
        preselected_kmer_counts = Counter()

    # process each graph component separately
    component_iter = tqdm(
        nx.weakly_connected_components(dbg),
        unit="comp",
        desc="dbg components",
        total=num_components,
    )
    for component in component_iter:
        component_orfs = setreduce_attr(dbg, component, "orf")

        # generate all candidate tiles
        tile_to_name = {}
        for name in component_orfs:
            # special case short orfs
            if len(orfs[name]) < tile_size:
                tile_to_name.setdefault(orfs[name], []).append(
                    (name, 0, len(orfs[name]))
                )
            for (i, j, tile) in tile_op(orfs[name], tile_size, tile_size - 1):
                tile_to_name.setdefault(tile, []).append((name, i, j))
        candidate_tiles = list(tile_to_name.keys())
        tile_to_idx = dict([(tile, i) for (i, tile) in enumerate(candidate_tiles)])

        # generate init tile scores
        tile_scores = []
        tile_lens = []
        kmer_to_idxs = {}
        for idx, tile in enumerate(candidate_tiles):
            score = 0
            for kmer in set(gen_kmers(tile, kmer_size)):
                score += dbg.nodes[kmer]["multiplicity"]
                kmer_to_idxs.setdefault(kmer, set()).add(idx)
            tile_scores.append(score / len(tile))
            tile_lens.append(len(tile))
        tile_scores = np.ma.asarray(tile_scores)
        tile_scores.harden_mask()
        tile_lens = np.asarray(tile_lens)

        # update tile scores with previously selected tiles and mask already
        # chosen tiles
        for kmer in set(preselected_kmer_counts.keys()) & set(kmer_to_idxs.keys()):
            idxs = list(kmer_to_idxs[kmer])
            tile_scores.data[idxs] -= (
                preselected_kmer_counts[kmer] * dbg.nodes[kmer]["multiplicity"]
            ) / len(tile)
        for tile in set(preselected_tiles) & set(candidate_tiles):
            tile_scores[tile_to_idx[tile]] = np.ma.masked

        # set upper bound on number of tiles for this component
        if kmer_cov:
            num_component_tiles = ceil(
                len(component) * kmer_cov / (tile_size - kmer_size + 1)
            )
        if num_tiles:
            num_component_tiles = ceil(
                len(component) / kmers_remaining * tiles_remaining
            )
            kmers_remaining -= len(component)

        # choose tiles
        for _ in range(num_component_tiles):
            # break out early if we've already sampled every tile
            if tile_scores.mask.sum() == len(tile_scores):
                break
            idx = tile_scores.argmax()
            tile_scores[idx] = np.ma.masked
            tile = candidate_tiles[idx]

            # write tile
            name, i, j = tile_to_name[tile][0]
            nterm = (
                "|NTERM" if dbg.nodes[tile[:kmer_size]].get("start_node", False) else ""
            )
            cterm = (
                "|CTERM" if dbg.nodes[tile[-kmer_size:]].get("end_node", False) else ""
            )
            print(f">{name}|{i}-{j}{nterm}{cterm}\n{tile}", file=output)

            # update tile scores
            for kmer in set(gen_kmers(tile, kmer_size)):
                idxs = list(kmer_to_idxs[kmer])
                tile_scores.data[idxs] -= (
                    dbg.nodes[kmer]["multiplicity"] / tile_lens[idxs]
                )
            tiles_remaining -= 1


@cli.command(short_help="compute stats on de Bruijn graph")
@option(
    "-d",
    "--dbg-path",
    type=Path(exists=True, dir_okay=False),
    required=True,
    help="pickled dbg",
)
@option(
    "-p",
    "--tiles",
    type=Path(exists=True, dir_okay=False),
    required=True,
    help="input protein tiles fasta",
)
@option(
    "-r",
    "--orfs",
    type=Path(exists=True, dir_okay=False),
    required=True,
    help="input ORFs fasta",
)
@option(
    "-o",
    "--output",
    type=File("w"),
    default=sys.stdout,
    help="output YAML path (default stdout)",
)
def dbgtilesummary(dbg_path, tiles, orfs, output):
    """Compute summary statistics on a de Bruin graph (DBG)

    The output is written out as a YAML. This operation computes statistics of
    an ORF set and a tile set relative to a DBG, which is why it requires those
    files too.

    """
    try:
        import networkx as nx
        from pepsyn.dbg import dbg_stats
    except ImportError:
        print("dbgtilesummary requires NetworkX", file=sys.stderr)
        raise Abort()

    dbg = nx.read_gpickle(dbg_path)
    with open(orfs, "r") as ip:
        orfs = [seq for (name, seq, qual) in readfq(ip)]
    with open(tiles, "r") as ip:
        tiles = [seq for (name, seq, qual) in readfq(ip)]

    stats = dbg_stats(dbg, orfs, tiles)

    print(yaml.dump(stats), file=output)


@cli.command(short_help="compute stats on peptides")
@option(
    "-p",
    "--tiles",
    type=Path(exists=True, dir_okay=False),
    required=True,
    help="input protein tiles fasta",
)
@option(
    "-r",
    "--orfs",
    type=Path(exists=True, dir_okay=False),
    required=True,
    help="input ORFs fasta",
)
@option(
    "-o",
    "--output",
    type=File("w"),
    default=sys.stdout,
    help="output YAML path (default stdout)",
)
def tilesummary(tiles, orfs, output):
    """Compute summary statistics on a set of peptide tiles

    These statistics are computed relative to a set of ORFs. This operation can
    be used on raw or cleaned ORFs.

    """
    with open(orfs, "r") as ip:
        orfs = {name: seq for (name, seq, qual) in readfq(ip)}
    with open(tiles, "r") as ip:
        tiles = {name: seq for (name, seq, qual) in readfq(ip)}
    stats = tile_stats(orfs, tiles)
    print(yaml.dump(stats), file=output)


@cli.command(short_help="compute stats on ORFs")
@argument_input
@argument_output
def orfsummary(input, output):
    """Compute summary statistics on a set of ORFs

    Can be used on raw or cleaned orfs.

    INPUT and OUTPUT are paths to fasta files or "-" to specify STDIN/STDOUT.

    """
    orfs = {name: seq for (name, seq, qual) in readfq(input)}
    stats = orf_stats(orfs)
    print(yaml.dump(stats), file=output)
