# config['library_name'] =
# config['input_orfs_path'] =
# config['tile_size'] =
# config['kmer_size'] =
# config['num_tiles'] =


rule all:
    input:
        'peptide_tiles.fasta',


rule cterm_pep_with_stops:
    input:
        config['input_orfs_path']
    output:
        'cterm_tiles.fasta'
    params:
        tile_size = config['tile_size']
    shell:
        """
        pepsyn --version
        cat {input} \
            | pepsyn ctermpep -l {params.tile_size} --add-stop - - \
            > {output}
        """


rule cluster_cterm_pep:
    input:
        'cterm_tiles.fasta'
    output:
        'cterm_tiles.clustered.fasta',
        'cterm_tiles.clustered.fasta.clstr'
    shell:
        """
        cd-hit -i {input} -o {output[0]} -c 0.95 -G 0 -aL 1.0 -aS 1.0 -M 0 -T 1 -d 0
        """


rule build_dbg:
    input:
        config['input_orfs_path']
    output:
        temp('dbgs/dbg.k{kmer_size}.pickle')
    params:
        kmer_size = '{kmer_size}'
    shell:
        """
        pepsyn --version
        pepsyn builddbg -k {params.kmer_size} {input} {output}
        """


rule dbg_sample_tiles:
    input:
        config['input_orfs_path'],
        'cterm_tiles.clustered.fasta',
        'dbgs/dbg.k{kmer_size}.pickle'
    output:
        'dbg_tiles.k{kmer_size}.fasta'
    params:
        tile_size = config['tile_size'],
        num_tiles = config['num_tiles']
    shell:
        """
        pepsyn --version
        NUM_CTERM=$(grep "^>" {input[1]} | wc -l)
        REMAINING_TILES=$(expr {params.num_tiles} - $NUM_CTERM)
        cat {input[0]} \
            | pepsyn greedykmercov -t {params.tile_size} -d {input[2]} \
                -n $REMAINING_TILES -p {input[1]} - - \
            > {output}
        """


rule dbg_pad_concat_tiles:
    input:
        expand('dbg_tiles.k{kmer_size}.fasta', kmer_size=config['kmer_size']),
        'cterm_tiles.clustered.fasta'
    output:
        'peptide_tiles.fasta'
    params:
        tile_size = config['tile_size']
    shell:
        """
        pepsyn --version
        cat {input} \
            | pepsyn pad -l {params.tile_size} --c-term - - \
            > {output}
        """
