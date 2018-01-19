# config['library_name'] =
# config['input_orfs_path'] =
# config['peptide_tiles_path'] =


ks = [5, 7, 9, 12, 15, 18, 21, 25]


rule all:
    input:
        expand('stats/dbg_stats.k{kmer_size}.yaml', kmer_size=ks),
        'report.html'


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


rule orf_stats:
    input:
        config['input_orfs_path']
    output:
        'stats/orf_stats.yaml'
    shell:
        """
        pepsyn --version
        pepsyn orfsummary {input} {output}
        """


rule tile_stats:
    input:
        config['peptide_tiles_path'],
        config['input_orfs_path']
    output:
        'stats/tile_stats.yaml'
    shell:
        """
        pepsyn --version
        pepsyn tilesummary -p {input[0]} -r {input[1]} -o {output}
        """


rule dbg_stats:
    input:
        'dbgs/dbg.k{kmer_size}.pickle',
        config['input_orfs_path'],
        config['peptide_tiles_path']
    output:
        'stats/dbg_stats.k{kmer_size}.yaml'
    shell:
        """
        pepsyn --version
        pepsyn dbgtilesummary -d {input[0]} -r {input[1]} -p {input[2]} -o {output}
        """


rule dbg_report_html:
    input:
        orf_stats = 'stats/orf_stats.yaml',
        tile_stats = 'stats/tile_stats.yaml',
        dbg_stats = expand('stats/dbg_stats.k{kmer_size}.yaml', kmer_size=ks)
    output:
        'report.html'
    run:
        from textwrap import dedent
        import yaml
        from bokeh.plotting import figure
        from bokeh.embed import components
        from bokeh.layouts import row
        from bokeh.resources import CDN as bokeh_CDN

        # load stats
        with open(input.orf_stats, 'r') as ip:
            orf_stats = yaml.load(ip)
        with open(input.tile_stats, 'r') as ip:
            tile_stats = yaml.load(ip)
        dbg_stats = []
        for stats_file in input.dbg_stats:
            with open(stats_file, 'r') as ip:
                dbg_stats.append(yaml.load(ip))
        dbg_stats.sort(key=lambda s: s['kmer_size'])
        k15_stats = next(filter(lambda s: s['kmer_size'] == 15, dbg_stats))

        # utility fns
        def extract_vals(attr):
            return [s[attr] for s in dbg_stats]

        def table_string(attr, fmt=None):
            fmt_str = '<td>{{:{}}}</td>'.format(fmt) if fmt else '<td>{}</td>'
            return ''.join([fmt_str.format(val) for val in extract_vals(attr)])

        def topk_string(attr):
            return '\n'.join(['            <tr><td><code>{}</code></td><td>{:,.1f}</td></tr>'.format(n, v) for n, v in tile_stats[attr]])

        # generate plots
        d = tile_stats['orf_coverage_hist']
        p1 = figure(plot_height=300, plot_width=300, title='ORF Coverage', x_axis_label='coverage')
        p1.quad(d['bin_edges'][:-1], d['bin_edges'][1:], d['hist'], 0)
        p2 = figure(plot_height=300, plot_width=300, title='ORF Coverage (log scale)', x_axis_label='coverage', y_axis_type='log')
        l, r, t = list(zip(*filter(lambda t: t[2] > 0, zip(d['bin_edges'][:-1], d['bin_edges'][1:], d['hist']))))
        p2.quad(l, r, t, 0.8)
        orf_coverage_hist = '\n'.join(components(row(p1, p2)))

        d = orf_stats['orf_lens_hist']
        p = figure(plot_height=300, plot_width=300, title='ORF lengths', x_axis_label='length (aa)')
        p.quad(d['bin_edges'][:-1], d['bin_edges'][1:], d['hist'], 0)
        orf_len_hist = '\n'.join(components(p))

        d = orf_stats['ambig_factor_hist']
        p = figure(plot_height=300, plot_width=300, title='Raw ORF ambiguities', x_axis_label='Number ORFS when disambiguated')
        p.quad(d['bin_edges'][:-1], d['bin_edges'][1:], d['hist'], 0)
        orf_ambig_factor_hist = '\n'.join(components(p))

        d = k15_stats['mult_hist']
        p1 = figure(plot_height=300, plot_width=300, title='k-mer multiplicity', x_axis_label='k-mer multiplicity')
        p1.quad(d['bin_edges'][:-1], d['bin_edges'][1:], d['hist'], 0)
        p2 = figure(plot_height=300, plot_width=300, title='k-mer multiplicity (log scale)', x_axis_label='k-mer multiplicity', y_axis_type='log')
        # filter out zeros for log scale
        l, r, t = list(zip(*filter(lambda t: t[2] > 0, zip(d['bin_edges'][:-1], d['bin_edges'][1:], d['hist']))))
        p2.quad(l, r, t, 0.8)
        mult_hist='\n'.join(components(row(p1, p2)))

        d = k15_stats['cov_hist']
        p1 = figure(plot_height=300, plot_width=300, title='k-mer coverage', x_axis_label='k-mer coverage')
        p1.quad(d['bin_edges'][:-1], d['bin_edges'][1:], d['hist'], 0)
        p2 = figure(plot_height=300, plot_width=300, title='k-mer coverage (log scale)', x_axis_label='k-mer coverage', y_axis_type='log')
        # filter out zeros for log scale
        l, r, t = list(zip(*filter(lambda t: t[2] > 0, zip(d['bin_edges'][:-1], d['bin_edges'][1:], d['hist']))))
        p2.quad(l, r, t, 0.8)
        cov_hist='\n'.join(components(row(p1, p2)))

        rendered = dedent(
            """
            <!DOCTYPE html>
            <html lang="en">
                <head>
                    <meta charset="utf-8">
                    <link rel="stylesheet" href="{bokeh_css_url}" type="text/css" />
                    <script type="text/javascript" src="{bokeh_js_url}"></script>
                    <title>{library_name}</title>
                    <style>
                        html {{font-family: 'helvetica neue', helvetica, arial, sans-serif;}}
                        body {{margin: 50px;}}
                        table {{border-collapse: collapse; border-top: 1px solid; border-bottom: 1px solid; margin-top: 20px; margin-bottom: 20px;}}
                        thead {{border-collapse: collapse; border-bottom: 1px solid;}}
                        th, td {{padding: 5px;}}
                        tr td {{text-align: right;}}
                        tr td:nth-child(1) {{text-align: left; width: 180px; padding-left: 20px; text-indent: -20px}}
                        .bk-root {{margin-top: 10px; margin-bottom: 20px;}}
                    </style>
                </head>
                <body>
                    <h1>{library_name}</h1>

                    <h2>Design summary</h2>

                    <table>
                        <tr><td>Tile size</td><td>{tile_size}</td></tr>
                        <tr><td>Number of tiles</td><td>{num_tiles:,}</td></tr>
                        <tr><td>Total tile residues</td><td>{total_tile_residues:,d}</td></tr>
                        <tr><td>Average ORF coverage</td><td>{avg_orf_coverage:,.1f}×</td></tr>
                        <tr><td>Num ORFs smaller than tile size ({tile_size})</td><td>{num_orfs_smaller_than_tile_size:,d}</td></tr>
                    </table>

                    {orf_coverage_hist}

                    <h2>Input ORFs</h2>

                    <table>
                        <tr><td>Number of raw input ORFs</td><td>{num_raw_orfs:,d}</td></tr>
                        <tr><td>Number of cleaned ORFs</td><td>{num_clean_orfs:,d}</td></tr>
                        <tr><td>Total ORF residues (disambig)</td><td>{total_orf_residues:,d}</td></tr>
                        <tr><td>Max ambiguity factor</td><td>{max_ambig_factor:,.1f}</td></tr>
                    </table>

                    {orf_len_hist}

                    {orf_ambig_factor_hist}

                    <table>
                        <thead>
                            <tr><th>Most-covered ORFs</th><th>coverage</th></tr>
                        </thead>
                        {top_5_orf_cov}
                    </table>

                    <table>
                        <thead>
                            <tr><th>Least-covered ORFs</th><th>coverage</th></tr>
                        </thead>
                        {bot_5_orf_cov}
                    </table>


                    <h2>De Bruijn graph (DBG) analysis</h2>

                    <p>Properties of the DBGs for different values of k.</p>

                    <table>
                        <thead>
                            <tr><th>k</th>{ks}</tr>
                        </thead>
                        <tr><td>Number of observed k-mers</td>{num_observed_kmers}</tr>
                        <tr><td>Num mult-1 kmers</td>{num_multiplicity_1_kmers}</tr>
                        <tr><td>Num mult > 1 kmers</td>{num_multiplicity_gt1_kmers}</tr>
                        <tr><td>Num k-mers covered</td>{num_kmers_covered}</tr>
                        <tr><td>Max k-mer multiplicity</td>{max_kmer_multiplicity}</tr>
                        <tr><td>Avg k-mer multiplicity</td>{avg_kmer_multiplicity}</tr>
                        <tr><td>Avg k-mer coverage</td>{avg_kmer_coverage}</tr>
                        <tr><td>Median k-mer coverage</td>{median_kmer_coverage}</tr>
                        <tr><td>Max k-mer coverage</td>{max_kmer_coverage}</tr>
                        <tr><td>Num k-mers missed</td>{num_kmers_missed}</tr>
                        <tr><td>Frac k-mers covered</td>{frac_kmers_covered}</tr>
                        <tr><td>Frac k-mers coverage > 1</td>{frac_kmers_covered_gt1}</tr>
                        <tr><td>C-term k-mer coverage</td>{cterm_kmer_cov}</tr>
                        <tr><td>N-term k-mer coverage</td>{nterm_kmer_cov}</tr>
                        <tr><td>Frac of mult-weighted kmers covered</td>{frac_mult_weighted_kmers_covered}</tr>
                        <tr><td>Number component subgraphs in DBG</td>{num_dbg_components}</tr>
                    </table>

                    <h3>Multiplicity and coverage histograms for 15-mer DBG</h3>

                    {mult_hist}
                    {cov_hist}
                </body>
            </html>
            """).format(
                bokeh_css_url=bokeh_CDN.css_files[0],
                bokeh_js_url=bokeh_CDN.js_files[0],
                library_name=config['library_name'],
                tile_size=tile_stats['tile_size'],
                num_tiles=tile_stats['num_tiles'],
                total_tile_residues=tile_stats['total_tile_residues'],
                avg_orf_coverage=tile_stats['avg_orf_coverage'],
                num_orfs_smaller_than_tile_size=tile_stats['num_orfs_smaller_than_tile_size'],
                orf_coverage_hist=orf_coverage_hist,
                num_raw_orfs=orf_stats['num_orfs'],
                num_clean_orfs=orf_stats['num_orfs'],
                total_orf_residues=orf_stats['total_orf_residues'],
                max_ambig_factor=orf_stats['max_ambig_factor'],
                orf_len_hist=orf_len_hist,
                orf_ambig_factor_hist=orf_ambig_factor_hist,
                top_5_orf_cov=topk_string('top_5_orf_cov'),
                bot_5_orf_cov=topk_string('bot_5_orf_cov'),
                ks=''.join(['<th>{}</th>'.format(v) for v in extract_vals('kmer_size')]),
                num_observed_kmers=table_string('num_observed_kmers', ',d'),
                max_kmer_multiplicity=table_string('max_kmer_multiplicity', ',d'),
                avg_kmer_multiplicity=table_string('avg_kmer_multiplicity', ',.2f'),
                num_multiplicity_1_kmers=table_string('num_multiplicity_1_kmers', ',d'),
                num_multiplicity_gt1_kmers=table_string('num_multiplicity_gt1_kmers', ',d'),
                avg_kmer_coverage=table_string('avg_kmer_coverage', ',.2f'),
                median_kmer_coverage=table_string('median_kmer_coverage', ',d'),
                max_kmer_coverage=table_string('max_kmer_coverage', ',d'),
                num_kmers_covered=table_string('num_kmers_covered', ',d'),
                num_kmers_missed=table_string('num_kmers_missed', ',d'),
                frac_kmers_covered=table_string('frac_kmers_covered', '.2f'),
                frac_kmers_covered_gt1=table_string('frac_kmers_covered_gt1', '.2f'),
                cterm_kmer_cov=table_string('cterm_kmer_cov', '.2f'),
                nterm_kmer_cov=table_string('nterm_kmer_cov', '.2f'),
                frac_mult_weighted_kmers_covered=table_string('frac_mult_weighted_kmers_covered', '.2f'),
                num_dbg_components=table_string('num_dbg_components', ',d'),
                mult_hist=mult_hist,
                cov_hist=cov_hist)
        with open(output[0], 'w') as op:
            print(rendered, file=op)


rule rule dbg_report_pdf:
    input:
        'report.html'
    output:
        'report.pdf'
    shell:
        """
        /Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome \
            --headless --virtual-time-budget=1000 --disable-gpu --print-to-pdf={output} {input}
        """
