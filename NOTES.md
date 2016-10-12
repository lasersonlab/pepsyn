**Notes from Church lab folks (DBG, Gleb, Marc Lajoie)**

* Start with wild type sequence and mutate N-terminal codons to remove secondary structure.
* Gleb has Python modules that overfit on their Recoli-57 project
* https://github.com/churchlab/getk/blob/master/src/getk/codon_usage_memex.py
    * Factory caller at the bottom `get_ecoli_codon_usage_memex()` will give you back a usage memex object.
* A messier version but with more codon usage calc
    * https://github.com/churchlab/genome-recoder/blob/master/src/codon_usage_memex.py
* Since you’re probably aiming for high expression, you probably just want to write a function to predict min free energy of your 5-prime mRNA
    * You can do so using https://github.com/churchlab/getk/blob/master/src/getk/hybrid_ss_min_wrapper.py, which requires installing Unafold and specifically the hybrid_ss_min tool
    * We typically pass (-30, +100) relative to start codon
    * So generate some sequences and pick one for which predicted MFE is closer to 0 (less negative) as possible
    * NOTE the https://github.com/churchlab/getk repo is still private and we have a license/patent situation going on so let me know if you pull out any parts into a separate open source library.
* DBG and Gleb still think it would be beneficial to get a real nucleotide sequence and tweak that rather than generate a sequence from scratch.  Even if you have to tweak it a lot to fix GC content / codon usage.
* DBG: If you are moving from low GC to high GC, then wholesale codon opt makes more sense, but otherwise, I would just keep it WT and optimize the n-term folding.
* Marc: as a more general answer for arbitrary genes, you can also use dnaworks (I think that's open source and works on command line).
