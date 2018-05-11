from glob import glob

from tqdm import tqdm

from llutil.utils import fastx_stem


# CONFIGURE
config['library_name'] = 'staph'
config['input_glob'] = '/Users/laserson/Downloads/BaseSpace/phip-15-65368304/FASTQ_Generation_2018-02-15_07_20_21Z-79385102/T7Staph_L001-ds.e545beacf4ca40a9b2e84547e9cabf0c/*.fastq.gz'
config['reference_index'] = '/Users/laserson/ll-dropbox/phage_libraries_private/staph/indexes/bowtie2/staph'


input_files = glob(config['input_glob'])
config['samples'] = {fastx_stem(f): f for f in input_files}


rule all:
    input:
        expand('{sample}.bam', sample=config['samples']),
        expand('{sample}.counts.tsv', sample=config['samples']),
        expand('{sample}.summary.yaml', sample=config['samples']),
        expand('{sample}.report.html', sample=config['samples']),


rule align_ref:
    input:
        fastq = lambda wildcards: config['samples'][wildcards.sample]
    output:
        '{sample}.bam',
    params:
        reference_index = config['reference_index']
    threads: 4
    shell:
        """
        bowtie2 -p {threads} -x {params.reference_index} -U {input} \
                   | samtools sort -O BAM -o {output}
        """

rule compute_counts:
    input:
        '{sample}.bam'
    output:
        '{sample}.counts.tsv'
    shell:
        """
        echo "id\tcount" > {output}
        samtools depth -aa -m 1000000 {input} \
            | tqdm \
            | awk 'BEGIN {{OFS="\t"}} {{counts[$1] = ($3 < counts[$1]) ? counts[$1] : $3}} END {{for (c in counts) {{print c, counts[c]}}}}' \
            >> {output}
        """


rule compute_snps:
    input:
        '{sample}.bam'
    output:
        '{sample}.varscan_snps.tsv'
    params:
        reference_index = config['reference_index']
    shell:
        """
        samtools mpileup -f {params.reference_index} {input} \
            | varscan mpileup2snp \
            > {output}
        """


rule compute_indels:
    input:
        '{sample}.bam'
    output:
        '{sample}.varscan_indels.tsv'
    params:
        reference_index = config['reference_index']
    shell:
        """
        samtools mpileup -f {params.reference_index} {input} \
            | varscan mpileup2indel \
            > {output}
        """


rule compute_summary:
    input:
        '{sample}.bam'
    output:
        '{sample}.summary.yaml'
    run:
        from pysam import AlignmentFile
        import yaml

        s = {}
        s['num_reads'] = 0
        s['num_reads_aligned'] = 0
        s['num_reads_aligned_perfect'] = 0
        clones_aligned = set()
        clones_aligned_perfect = set()

        bamfile = AlignmentFile(input[0], 'rb')
        for aln in tqdm(bamfile):
            s['num_reads'] += 1
            if not aln.is_unmapped:
                s['num_reads_aligned'] += 1
                clones_aligned.add(aln.reference_name)
                if len(aln.cigar) == 1 and aln.cigar[0] == (0, aln.query_length):
                    s['num_reads_aligned_perfect'] += 1
                    clones_aligned_perfect.add(aln.reference_name)

        s['read_length'] = aln.query_length
        s['num_clones_with_alignment'] = len(clones_aligned)
        s['num_clones_with_alignment_perfect'] = len(clones_aligned_perfect)

        with open(output[0], 'w') as op:
            print(yaml.dump(s), file=op)


rule generate_report:
    input:
        stats = '{sample}.summary.yaml',
        counts = '{sample}.counts.tsv',
    output:
        '{sample}.report.html'
    run:
        from textwrap import dedent

        import yaml
        import numpy as np
        import pandas as pd
        from bokeh.plotting import figure
        from bokeh.embed import components
        from bokeh.layouts import row
        from bokeh.resources import CDN as bokeh_CDN

        s = {}
        s['library_name'] = config['library_name']
        s['bokeh_css_url'] = bokeh_CDN.css_files[0]
        s['bokeh_js_url'] = bokeh_CDN.js_files[0]

        with open(input.stats, 'r') as ip:
            s.update(yaml.load(ip))
        counts = pd.read_csv(input.counts, sep='\t', header=0)

        s['num_clones'] = len(counts)
        s['average_coverage'] = s['num_reads'] / s['num_clones']

        x0, x1, y0, y1 = list(range(len(counts))), list(range(len(counts))), 0, sorted(counts['count'])
        p1 = figure(plot_height=400, plot_width=400, title='clone sizes', x_axis_label='clones', y_axis_label='read count')
        p1.segment(x0, y0, x1, y1)
        # filter out zeros for log scale
        x0nz, x1nz, y1nz = list(zip(*filter(lambda t: t[2] > 0, zip(x0, x1, y1))))
        p2 = figure(plot_height=400, plot_width=400, title='clone sizes (log scale)', x_axis_label='clones', y_axis_label='read count', y_axis_type='log')
        p2.segment(x0nz, 0.1, x1nz, y1nz)
        p2.x_range.start = 0
        s['clone_sizes'] = '\n'.join(components(row(p1, p2)))

        hist, edges = np.histogram(counts['count'], bins=range(counts['count'].max()))
        p3 = figure(plot_height=400, plot_width=400, title='coverage', x_axis_label='read count', y_axis_label='num clones', x_range=(0, 100), y_range=(0, 1000))
        p3.quad(edges[:-1], edges[1:], hist, 0)
        # filter out zeros for log scale
        l, r, t = list(zip(*filter(lambda t: t[2] > 0, zip(edges[:-1], edges[1:], hist))))
        p4 = figure(plot_height=400, plot_width=400, title='coverage (log scale)', x_axis_label='read count', y_axis_label='num clones', y_axis_type='log')
        p4.quad(l, r, t, 0.8)
        s['count_hist'] = '\n'.join(components(row(p3, p4)))

        hist, edges = np.histogram(counts['count'], bins=range(counts['count'].max()), density=True)
        cumsum = [0.] + list(np.cumsum(hist))
        p5 = figure(plot_height=400, plot_width=400, title='coverage cdf', x_axis_label='read count', y_axis_label='cdf')
        p5.step(edges, cumsum, mode='after')
        p6 = figure(plot_height=400, plot_width=400, title='coverage cdf (log scale)', x_axis_label='read count', y_axis_label='cdf', y_axis_type='log')
        p6.step(edges[1:], cumsum[1:], mode='after')
        s['count_cdf'] = '\n'.join(components(row(p5, p6)))

        topk_freq_clones = [(tup.id, tup.count) for tup in counts.sort_values('count', ascending=False).iloc[:20].itertuples()]
        s['topk_string'] = '\n'.join(['<tr><td><code>{}</code></td><td>{:,d}</td></tr>'.format(x, y) for x, y in topk_freq_clones])

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
                        tr td:nth-child(1) {{text-align: left; width: 300px; padding-left: 20px; text-indent: -20px}}
                        .bk-root {{margin-top: 10px; margin-bottom: 20px;}}
                    </style>
                </head>
                <body>
                    <h1><code>{library_name}</code> synthetic library sequencing</h1>

                    <table>
                        <tr><td>Read length</td><td>{read_length}</td></tr>
                        <tr><td>Number of reads</td><td>{num_reads:,}</td></tr>
                        <tr><td>Average coverage (reads per clone)</td><td>{average_coverage:.1f}×</td></tr>
                        <tr><td>Number of reads aligned</td><td>{num_reads_aligned:,d}</td></tr>
                        <tr><td>Number of reads aligned perfectly</td><td>{num_reads_aligned_perfect:,d}</td></tr>
                        <tr><td>Number of clones (designed)</td><td>{num_clones:,d}</td></tr>
                        <tr><td>Number of clones with alignment</td><td>{num_clones_with_alignment:,d}</td></tr>
                        <tr><td>Number of clones with perfect alignment</td><td>{num_clones_with_alignment_perfect:,d}</td></tr>
                    </table>

                    {clone_sizes}
                    {count_hist}
                    {count_cdf}

                    <table>
                        <thead>
                            <tr><th>Most-frequent clones</th><th>read count</th></tr>
                        </thead>
                        {topk_string}
                    </table>
                </body>
            </html>
            """).format(**s)
        with open(output[0], 'w') as op:
            print(rendered, file=op)





