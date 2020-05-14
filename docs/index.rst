``pepsyn``
==========

``pepsyn`` is a Python package for designing oligonucleotide-encoded peptide
libraries. It provides a command line interface and a Python API for performing
common operations on DNA or protein sequences. In addition to sequence
manipulation, ``pepsyn`` can build a de Bruijn graph for use in peptide tile
selection and can also generate reports and statistical summaries.

From the command line, a typical pepsyn design "protocol" might look like:

.. code-block:: shell

    cat tests/proteins.fasta \
        | pepsyn tile -l 10 -p 3 - - \
        | pepsyn revtrans - - \
        | pepsyn prefix -p ACGGG - - \
        | pepsyn suffix -s TGCTG - - \
        | pepsyn recodesite --site EcoRI --clip-left 5 --clip-right 5 - -

Tools can be chained together with pipes like standard POSIX tools. In the
example above, each protein sequence in ``proteins.fasta`` will be chopped up
into 10 amino acid tiles that overlap each other by 3 amino acids. Each of these
tiles is then reverse-translated into a DNA sequence, prefixed and suffixed with
adaptor sequences, and finally the coding region is recoded to ensure there are
no EcoRI restriction enzyme sites.

You can get the library directly from PyPI:

.. code-block:: shell

    pip install pepsyn

The only heavyweight required dependency is `Biopython <https://biopython.org/>`_.
To make use of the de Bruijn graph functionality, you additionally need
`NetworkX <https://networkx.github.io/>`_.

Documentation
-------------

.. toctree::
    :maxdepth: 2

    tutorial
    cli
    api
