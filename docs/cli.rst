``pepsyn`` command line reference
=================================


Overview
--------

.. click:: pepsyn.cli:cli
   :prog: pepsyn


Operations on any sequences
---------------------------

.. click:: pepsyn.cli:clip
    :prog: pepsyn clip

.. click:: pepsyn.cli:filterlen
    :prog: pepsyn filterlen

.. click:: pepsyn.cli:tile
    :prog: pepsyn tile



Operations on DNA sequences
---------------------------

.. click:: pepsyn.cli:findsite
    :prog: pepsyn findsite

.. click:: pepsyn.cli:prefix
    :prog: pepsyn prefix

.. click:: pepsyn.cli:suffix
    :prog: pepsyn suffix

.. click:: pepsyn.cli:recodesite
    :prog: pepsyn recodesite

.. click:: pepsyn.cli:translate
    :prog: pepsyn translate


Operations on protein sequences
-------------------------------

.. click:: pepsyn.cli:ctermpep
    :prog: pepsyn ctermpep

.. click:: pepsyn.cli:disambiguateaa
    :prog: pepsyn disambiguateaa

.. click:: pepsyn.cli:filterstop
    :prog: pepsyn filterstop

.. click:: pepsyn.cli:pad
    :prog: pepsyn pad

.. click:: pepsyn.cli:revtrans
    :prog: pepsyn revtrans

.. click:: pepsyn.cli:stripstop
    :prog: pepsyn stripstop

.. click:: pepsyn.cli:x2ggsg
    :prog: pepsyn x2ggsg



de Bruijn graph tools
^^^^^^^^^^^^^^^^^^^^^

.. click:: pepsyn.cli:builddbg
    :prog: pepsyn builddbg

.. click:: pepsyn.cli:greedykmercov
    :prog: pepsyn greedykmercov



Summary/report generation
-------------------------

.. click:: pepsyn.cli:orfsummary
    :prog: pepsyn orfsummary

.. click:: pepsyn.cli:tilesummary
    :prog: pepsyn tilesummary

.. click:: pepsyn.cli:dbgtilesummary
    :prog: pepsyn dbgtilesummary
