from glob import glob

from llutil.utils import fastx_stem


input_files = glob(config['input_glob'])
config['samples'] = {fastx_stem(f): f for f in input_files}


rule all:
    input:
        expand('species_aligned/{sample}.bam', sample=config['samples']),
        expand('annot/{sample}.annot_inserts.tsv', sample=config['samples']),
        expand('coverage/{sample}.cov_stranded.tsv', sample=config['samples']),
        expand('coverage/{sample}.cov_unstranded.tsv', sample=config['samples']),
        expand('coverage/{sample}.cov_inframe.tsv', sample=config['samples']),
        expand('coverage/{sample}.cov_genome.tsv', sample=config['samples']),
        expand('clustered_reads/{sample}.clustered.counts.tsv', sample=config['samples']),
        expand('annot/{sample}.inframe_inserts.counts.tsv', sample=config['samples']),
        # expand('diversity/{sample}.recon_allreads.txt', sample=config['samples']),
        # expand('diversity/{sample}.recon_inframe.txt', sample=config['samples'])


rule deplete_vector:
    input:
        fastq = lambda wildcards: config['samples'][wildcards.sample]
    output:
        fastq = 'vector_depleted/{sample}.fastq',
        tmpdir = temp('tmp_{sample}')
    params:
        insert_length = config['insert_length'],
        insert_stddev = config['insert_stddev'],
        vector_index = config['vector_kallisto_index']
    shell:
        """
        echo REENABLE MULTITHREADED KALLISTO
        echo REENABLE MULTITHREADED KALLISTO
        echo REENABLE MULTITHREADED KALLISTO
        echo REENABLE MULTITHREADED KALLISTO
        echo REENABLE MULTITHREADED KALLISTO

        kallisto pseudo --single --pseudobam -t 1 -l {params.insert_length} -s {params.insert_stddev} -o {output.tmpdir} -i {params.vector_index} {input} \
            | samtools fastq -f 4 - \
            > {output.fastq}
        """


rule align_species:
    input:
        'vector_depleted/{sample}.fastq'
    output:
        'species_aligned/{sample}.bam'
    params:
        species_index = config['species_index']
    threads: 4
    shell:
        """
        bowtie2 -p {threads} -x {params.species_index} -U {input} \
                   | samtools sort -@ {threads} -O BAM -o {output}
        """


rule gene_filter_gff:
    input:
        config['species_gff']
    output:
        temp('genes.gff')
    shell:
        r"""cat {input} | awk -F \t '/^[^#]/ && $3 == "gene"' > {output}"""


rule contig_lengths:
    input:
        config['species_fasta']
    output:
        temp('chrom_sizes.tsv')
    shell:
        """
        samtools dict -H {input} \
            | awk 'BEGIN {{FS="\t"; OFS="\t"}} {{print substr($2, 4), substr($3, 4)}}' \
            > {output}
        """


rule generate_inserts:
    input:
        reads = 'species_aligned/{sample}.bam',
        genome = 'chrom_sizes.tsv'
    output:
        temp('{sample}.inserts.bed')
    params:
        extend = config['insert_length'] - config['read_length']
    shell:
        """
        bedtools bamtobed -cigar -i {input.reads} \
            | bedtools slop -i stdin -s -g {input.genome} -l 0 -r {params.extend} \
            > {output}
        """


rule annotate_inserts_with_genes:
    input:
        inserts = rules.generate_inserts.output,
        genes = 'genes.gff'
    output:
        'annot/{sample}.annot_inserts.tsv'
    shell:
        'bedtools intersect -s -a {input.inserts} -b {input.genes} -wao > {output}'


rule filter_inframe_inserts:
    input:
        'annot/{sample}.annot_inserts.tsv'
    output:
        'annot/{sample}.inframe_inserts.tsv'
    shell:
        """
        # first awk script filters reads that overlap something
        # second awk script filters in-frame
        cat {input} \
            | awk '$11 >= 0' \
            | awk '($6 == "+" && ($2 - $11 + 1) % 3 == 0) || ($6 == "-" && ($3 - $12) % 3 == 0)' \
            > {output}
        """


rule stranded_gene_coverage:
    input:
        inserts = rules.generate_inserts.output,
        genes = 'genes.gff'
    output:
        'coverage/{sample}.cov_stranded.tsv'
    shell:
        'bedtools coverage -s -a {input.genes} -b {input.inserts} -wao > {output}'


rule unstranded_gene_coverage:
    input:
        inserts = rules.generate_inserts.output,
        genes = 'genes.gff'
    output:
        'coverage/{sample}.cov_unstranded.tsv'
    shell:
        'bedtools coverage -a {input.genes} -b {input.inserts} -wao > {output}'


rule genome_coverage:
    input:
        inserts = rules.generate_inserts.output,
        genome = 'chrom_sizes.tsv'
    output:
        'coverage/{sample}.cov_genome.tsv'
    shell:
        'bedtools genomecov -d -g {input.genome} -i {input.inserts} > {output}'


rule inframe_gene_coverage:
    input:
        inserts = 'annot/{sample}.inframe_inserts.tsv',
        genes = 'genes.gff'
    output:
        'coverage/{sample}.cov_inframe.tsv'
    shell:
        'bedtools coverage -s -a {input.genes} -b {input.inserts} -wao > {output}'



rule clustered_read_counts:
    input:
        'vector_depleted/{sample}.fastq'
    output:
        'clustered_reads/{sample}.clustered.fastq',
        'clustered_reads/{sample}.clustered.fastq.clstr',
        'clustered_reads/{sample}.clustered.counts.tsv'
    params:
        word_size = 8,
        threshold = 0.95
    threads: 4
    shell:
        """
        cd-hit-est -i {input} -o {output[0]} -n {params.word_size} -c {params.threshold} -r 0 -d 0 -M 0 -T {threads}
        llutil cd-hit collapse -r {output[0]} -c {output[1]} -o {output[2]}
        """


rule recon_diversity_allreads:
    input:
        'clustered_reads/{sample}.clustered.counts.tsv'
    output:
        temp('diversity/tmp_{sample}.recon_input.tsv'),
        'diversity/{sample}.recon_allreads.txt'
    params:
        python2 = config['python2'],
        recon = config['recon']
    shell:
        r"""
        cat {input} | awk '{{print $2 "\t" $1}}' > {output[0]}
        {params.python2} {params.recon} -R -t 1000 -o {output[1]} {output[0]}
        """


rule inframe_insert_counts:
    input:
        'annot/{sample}.inframe_inserts.tsv'
    output:
        'annot/{sample}.inframe_inserts.counts.tsv'
    shell:
        r"""
        # cut | sort | uniq -c generates the counts
        # the second sort | uniq filters out repeated reads, since they can overlap
        # multiple genes and generate multiple lines in the join
        cat {input} \
            | cut -f 1-3,8,11,12 | sort | uniq -c \
            | awk '{{print $2, $3, $4 "\t" $1}}' \
            | sort -t \t -k 1 | uniq \
            > {output}
        """


rule recon_diversity_inframe:
    input:
        'annot/{sample}.inframe_inserts.counts.tsv'
    output:
        'diversity/{sample}.recon_inframe.txt'
    params:
        python2 = config['python2'],
        recon = config['recon']
    shell:
        r"""
        {params.python2} {params.recon} -R -t 1000 -o {output} {input}
        """

