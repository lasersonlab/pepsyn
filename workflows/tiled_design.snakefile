# config['library_name'] =
# config['input_orfs_path'] =
# config['tile_size'] =
# config['overlap'] =


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


rule naive_tiling:
    input:
        config['input_orfs_path']
    output:
        'orf_tiles.fasta'
    params:
        tile_size = config['tile_size'],
        overlap = config['overlap']
    shell:
        """
        pepsyn --version
        cat {input} \
            | pepsyn tile -l {params.tile_size} -p {params.overlap} - - \
            > {output}
        """


rule cluster_naive_tiles:
    input:
        'orf_tiles.fasta'
    output:
        'orf_tiles.clustered.fasta',
        'orf_tiles.clustered.fasta.clstr'
    shell:
        """
        cd-hit -i {input} -o {output[0]} -c 0.95 -G 0 -A 50 -M 0 -T 1 -d 0
        """


rule naive_pad_concat_tiles:
    input:
        'orf_tiles.clustered.fasta',
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
