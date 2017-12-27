# `pepsyn`
Peptide library design


**Installation**

`pepsyn` is developed for Python 3.6+ and requires `biopython`.  The CLI uses
`click`.

```bash
pip install pepsyn
```

or for the latest development version

```bash
git clone https://github.com/lasersonlab/pepsyn.git
cd pepsyn
python setup.py install
```

To run the tests

```bash
py.test
```

**Usage**

`pepsyn` can be used as a library, but also comes with a CLI.

The main CLI command is `pepsyn` and it can read from `stdin` and `stdout` by
specifying `-` as the input or output.  Input and output files are FASTA
formatted.

```bash
$ pepsyn -h
Usage: pepsyn [OPTIONS] COMMAND [ARGS]...

  pepsyn -- peptide synthesis design

Options:
  -h, --help  Show this message and exit.

Commands:
  prefix      add a prefix to each sequence
  recodesite  remove site from each sequence's CDS by...
  revtrans    reverse translate amino acid sequences into...
  stats       compute some sequence statistics
  suffix      add a suffix to each sequence
  tile        tile a set of sequences

```

For example, to tile a set of sequences, use the `pepsyn tile` command.  Each
subcommand has its own help message with relevant options.

```bash
$ pepsyn tile -h
Usage: pepsyn tile [OPTIONS] INPUT OUTPUT

  tile a set of sequences

Options:
  -l, --length INTEGER   Length of output oligos
  -p, --overlap INTEGER  Overlap of oligos
  -h, --help             Show this message and exit.
```

Commands can be piped into each other.  (Note: `stdin` and `stdout` are
signified with `-`.)

```bash
cat pepsyn/tests/proteins.fasta \
    | pepsyn tile -l 10 -p 3 - - \
    | pepsyn revtrans - - \
    | pepsyn prefix -p ACGGG - - \
    | pepsyn suffix -s TGCTG - - \
    | pepsyn recodesite --site EcoRI --clip-left 5 --clip-right 5 - -
```
