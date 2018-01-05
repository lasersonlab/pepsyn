# config['library_name'] =
# config['input_orfs_path'] =
# config['tile_size'] =
# config['kmer_size'] =
# config['num_tiles'] =


ks = [5, 7, 9, 12, 15, 18, 21, 25]
if config['kmer_size'] not in ks:
    ks.append(config['kmer_size'])


rule all:
    input:
        'peptide_tiles.fasta',
        expand('stats/dbg_orf_tile_stats.k{kmer_size}.yaml', kmer_size=ks),
        'report.md',
        'report.pdf'


rule clean_input_orfs:
    input:
        config['input_orfs_path']
    output:
        'input_orfs.clean.fasta'
    shell:
        """
        pepsyn --version
        cat {input} \
            | pepsyn stripstop - - \
            | pepsyn filterstop - - \
            | pepsyn x2ggsg - - \
            | pepsyn disambiguateaa - - \
            > {output}
        """


rule cterm_pep_with_stops:
    input:
        'input_orfs.clean.fasta'
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
        cd-hit -i {input} -o {output[0]} -c 0.9 -G 0 -aL 1.0 -aS 1.0 -M 0 -T 3 -d 0
        """


rule build_dbg:
    input:
        'input_orfs.clean.fasta'
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
        'input_orfs.clean.fasta',
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


rule concat_tiles:
    input:
        expand('dbg_tiles.k{kmer_size}.fasta', kmer_size=config['kmer_size']),
        'cterm_tiles.clustered.fasta'
    output:
        'peptide_tiles.fasta'
    shell:
        """
        cat {input} > {output}
        """


rule dbg_orf_tile_stats:
    input:
        'dbgs/dbg.k{kmer_size}.pickle',
        'input_orfs.clean.fasta',
        'peptide_tiles.fasta'
    output:
        'stats/dbg_orf_tile_stats.k{kmer_size}.yaml'
    shell:
        """
        pepsyn --version
        pepsyn dbgtilesummary -d {input[0]} -r {input[1]} -p {input[2]} -o {output}
        """

rule dbg_report_md:
    input:
        expand('stats/dbg_orf_tile_stats.k{kmer_size}.yaml', kmer_size=ks)
    output:
        'figs/',
        'report.md'
    run:
        from textwrap import dedent

        import yaml
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        plt.style.use('seaborn')
        mpl.rcParams['figure.dpi'] = 300

        stats = []
        for stats_file in input:
            with open(stats_file, 'r') as ip:
                stats.append(yaml.load(ip))
        stats.sort(key=lambda s: s['kmer_size'])

        fig, ax = plt.subplots()
        bin_edges = stats[0]['cov_hist']['bin_edges']
        hist = stats[0]['cov_hist']['hist']
        ax.hlines(hist, bin_edges[:-1], bin_edges[1:])
        ax.set(xlim=(0, 100), xlabel='coverage', title='kmer coverage hist of tiling')
        fig.savefig('figs/cov_hist.png')

        fig, ax = plt.subplots()
        bin_edges = stats[0]['mult_hist']['bin_edges']
        hist = stats[0]['mult_hist']['hist']
        ax.hlines(hist, bin_edges[:-1], bin_edges[1:])
        ax.set(xlim=(0, 100), xlabel='multiplicity', title='kmer multiplicity hist of orfs')
        fig.savefig('figs/mult_hist.png')

        fig, ax = plt.subplots()
        bin_edges = stats[0]['orf_lens']['bin_edges']
        hist = stats[0]['orf_lens']['hist']
        ax.hlines(hist, bin_edges[:-1], bin_edges[1:])
        ax.set(xlabel='len', title='orf len hist')
        fig.savefig('figs/orf_lens.png')

        fig, ax = plt.subplots()
        bin_edges = stats[0]['tile_lens']['bin_edges']
        hist = stats[0]['tile_lens']['hist']
        ax.hlines(hist, bin_edges[:-1], bin_edges[1:])
        ax.set(xlabel='len', title='tile len hist', yscale='log')
        fig.savefig('figs/tile_lens.png')

        def extract_vals(attr):
            return [s[attr] for s in stats]

        def kmer_varying_line_fig(attr):
            fig, ax = plt.subplots()
            ax.plot(extract_vals('kmer_size'), extract_vals(attr))
            ax.set(xlabel='kmer size', ylabel=attr)
            fig.savefig(f'figs/{attr}.png')

        design_kmer_size_idx = extract_vals('kmer_size').index(config['kmer_size'])
        design_stats = stats[design_kmer_size_idx]

        kmer_varying_line_fig('avg_kmer_coverage')
        kmer_varying_line_fig('avg_kmer_multiplicity')
        kmer_varying_line_fig('cterm_kmer_cov')
        kmer_varying_line_fig('frac_kmers_covered')
        kmer_varying_line_fig('max_kmer_coverage')
        kmer_varying_line_fig('max_kmer_multiplicity')
        kmer_varying_line_fig('median_kmer_coverage')
        kmer_varying_line_fig('min_theor_tiles_1x_cov')
        kmer_varying_line_fig('frac_mult_weighted_kmers_covered')
        kmer_varying_line_fig('nterm_kmer_cov')
        kmer_varying_line_fig('num_dbg_components')
        kmer_varying_line_fig('num_kmers_covered')
        kmer_varying_line_fig('num_kmers_missed')
        kmer_varying_line_fig('num_multiplicity_1_kmers')
        kmer_varying_line_fig('num_observed_kmers')
        kmer_varying_line_fig('num_orfs_smaller_than_kmer_size')

        def table_string(attr, fmt=None):
            fmt_str = '{{:{}}}'.format(fmt) if fmt else '{}'
            return ' | '.join([fmt_str.format(val) for val in extract_vals(attr)])

        rendered = dedent(
            """
            # {library_name} - {kmer_size}-mer design


            ## Design summary

            |||
            |---|---|
            | k-mer size | {kmer_size} |
            | Tile size | {tile_size} |
            | Number of tiles | {num_tiles} |
            | Total tile residues | {total_tile_residues} |
            | Average ORF coverage | {avg_orf_coverage:.1f}× |
            | Approx num tiles naive 1x tiling | {approx_num_tiles_naive_1x_tiling} |
            | Max theor {kmer_size}-mers per tile | {max_theor_kmers_per_tile} |
            | Min theor tiles for 1x k-mer cov | {min_theor_tiles_1x_cov} |
            | Num ORFs smaller than k-mer size ({kmer_size}) | {num_orfs_smaller_than_kmer_size} |
            | Num ORFs smaller than tile size ({tile_size}) | {num_orfs_smaller_than_tile_size} |

            ![](figs/tile_lens.png)


            ## Input ORFs

            |||
            |---|---|
            | Number of ORFs | {num_orfs} |
            | Total ORF residues | {total_orf_residues} |

            ![](figs/orf_lens.png)


            ## De Bruijn graph (DBG) analysis

            Tile selection was based on the **{kmer_size}-mer** DBG.

            Properties of the DBGs for different values of `k`.

            | k | {ks} |
            |---|{header_sep}|
            | Number of observed k-mers | {num_observed_kmers} |
            | Max k-mer multiplicity | {max_kmer_multiplicity} |
            | Avg k-mer multiplicity | {avg_kmer_multiplicity} |
            | Num mult-1 kmers | {num_multiplicity_1_kmers} |
            | Num mult > 1 kmers | {num_multiplicity_gt1_kmers} |
            | Avg k-mer coverage | {avg_kmer_coverage} |
            | Median k-mer coverage | {median_kmer_coverage} |
            | Max k-mer coverage | {max_kmer_coverage} |
            | Num k-mers covered | {num_kmers_covered} |
            | Num k-mers missed | {num_kmers_missed} |
            | Frac k-mers covered | {frac_kmers_covered} |
            | Frac k-mers coverage > 1 | {frac_kmers_covered_gt1} |
            | C-term k-mer coverage | {cterm_kmer_cov} |
            | N-term k-mer coverage | {nterm_kmer_cov} |
            | Frac of mult-weighted kmers covered | {frac_mult_weighted_kmers_covered} |
            | Number component subgraphs in DBG | {num_dbg_components} |


            ### Multiplicity and coverage histograms for {kmer_size}-mer DBG

            ![](figs/mult_hist.png)

            ![](figs/cov_hist.png)
            """.format(
                header_sep='|'.join(['---'] * len(stats)),
                library_name=config['library_name'],
                kmer_size=config['kmer_size'],
                tile_size=config['tile_size'],
                num_tiles=design_stats['num_tiles'],
                total_tile_residues=design_stats['total_tile_residues'],
                avg_orf_coverage=design_stats['avg_orf_coverage'],
                approx_num_tiles_naive_1x_tiling=design_stats['approx_num_tiles_naive_1x_tiling'],
                max_theor_kmers_per_tile=design_stats['max_theor_kmers_per_tile'],
                min_theor_tiles_1x_cov=design_stats['min_theor_tiles_1x_cov'],
                num_orfs_smaller_than_kmer_size=design_stats['num_orfs_smaller_than_kmer_size'],
                num_orfs_smaller_than_tile_size=design_stats['num_orfs_smaller_than_tile_size'],
                num_orfs=design_stats['num_orfs'],
                total_orf_residues=design_stats['total_orf_residues'],
                ks=table_string('kmer_size'),
                num_observed_kmers=table_string('num_observed_kmers'),
                max_kmer_multiplicity=table_string('max_kmer_multiplicity'),
                avg_kmer_multiplicity=table_string('avg_kmer_multiplicity', '.2f'),
                num_multiplicity_1_kmers=table_string('num_multiplicity_1_kmers'),
                num_multiplicity_gt1_kmers=table_string('num_multiplicity_gt1_kmers'),
                avg_kmer_coverage=table_string('avg_kmer_coverage', '.2f'),
                median_kmer_coverage=table_string('median_kmer_coverage'),
                max_kmer_coverage=table_string('max_kmer_coverage'),
                num_kmers_covered=table_string('num_kmers_covered'),
                num_kmers_missed=table_string('num_kmers_missed'),
                frac_kmers_covered=table_string('frac_kmers_covered', '.3f'),
                frac_kmers_covered_gt1=table_string('frac_kmers_covered_gt1', '.3f'),
                cterm_kmer_cov=table_string('cterm_kmer_cov', '.2f'),
                nterm_kmer_cov=table_string('nterm_kmer_cov', '.2f'),
                frac_mult_weighted_kmers_covered=table_string('frac_mult_weighted_kmers_covered', '.2f'),
                num_dbg_components=table_string('num_dbg_components')))
        with open(output[1], 'w') as op:
            print(rendered, file=op)



rule rule dbg_report_pdf:
    input:
        'report.md'
    output:
        temp('report.html'),
        'report.pdf'
    shell:
        """
        pandoc {input} -f gfm -s -o {output[0]}
        pandoc {output[0]} -V geometry:margin=0.5in -o {output[1]}
        """
