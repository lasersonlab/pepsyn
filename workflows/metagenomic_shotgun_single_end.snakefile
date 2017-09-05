from glob import glob

from llutil.utils import fastx_stem


input_files = glob(config['input_glob'])
config['samples'] = {fastx_stem(f): f for f in input_files}

rule all:
    input:
        expand('centrifuge/{sample}.centrifuge_kreport.tsv', sample=config['samples']),
        expand('clustered_reads/{sample}.clustered.counts.tsv', sample=config['samples']),
        expand('clustered_reads/{sample}.aligned_clustered.counts.tsv', sample=config['samples']),


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


rule centrifuge_classify:
    input:
        'vector_depleted/{sample}.fastq'
    output:
        'centrifuge/{sample}.centrifuge_aln.tsv',
        'centrifuge/{sample}.centrifuge_kreport.tsv',
        temp('centrifuge/{sample}.centrifuge_report.tsv')
    params:
        index = config['centrifuge_index']
    threads: 4
    shell:
        """
        centrifuge -p {threads} -x {params.index} -U {input} -S {output[0]} --report-file {output[2]}
        centrifuge-kreport -x {params.index} {output[0]} > {output[1]}
        """

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


rule extract_aligned_reads:
    input:
        'centrifuge/{sample}.centrifuge_aln.tsv',
        'vector_depleted/{sample}.fastq',
    output:
        temp('tmp_aligned_reads.{sample}.fasta'),
    run:
        from Bio.SeqIO.QualityIO import FastqGeneralIterator
        ids = set()
        with open(input[0], 'r') as ip:
            _ = next(ip)
            for line in ip:
                fields = line.split('\t')
                if fields[1] != 'unclassified':
                    ids.add(line.split()[0])
        with open(input[1], 'r') as ip:
            with open(output[0], 'w') as op:
                for (title, seq, qual) in FastqGeneralIterator(ip):
                    if title in ids:
                        print('>{}\n{}'.format(title, seq), file=op)


rule clustered_aligned_read_counts:
    input:
        'tmp_aligned_reads.{sample}.fasta'
    output:
        'clustered_reads/{sample}.aligned_clustered.fasta',
        'clustered_reads/{sample}.aligned_clustered.fasta.clstr',
        'clustered_reads/{sample}.aligned_clustered.counts.tsv'
    params:
        word_size = 8,
        threshold = 0.95
    threads: 4
    shell:
        """
        cd-hit-est -i {input} -o {output[0]} -n {params.word_size} -c {params.threshold} -r 0 -d 0 -M 0 -T {threads}
        llutil cd-hit collapse -r {output[0]} -c {output[1]} -o {output[2]}
        """


# rule sketch_reads:
#     input:
#         'vector_depleted/{sample}.fastq'
#     output:
#         'sketches/{sample}.msh'
#     params:
#         min_kmer_copies = config['min_kmer_copies'],
#         size = 100000
#     threads: 4
#     shell:
#         'mash sketch -s {params.size} -p {threads} -m {params.min_kmer_copies} -o {output} {input}'


# rule mash_dist:
#     input:
#         'sketches/{sample}.msh'
#     output:
#         'dists/{sample}.dists.tsv'
#     params:
#         reference = config['reference_sketch']
#     threads: 4
#     shell:
#         'mash dist -p {threads} {params.reference} {input} > {output}'


# rule build_kallisto_index:
#     input:
#         'dists/{sample}.dists.tsv'
#     output:
#         'kallisto_idxs/{sample}.idx'
#     params:
#         top_n = config['top_n']
#     shell:
#         """
#         # SIGPIPE not handled properly with -o pipefail:
#         # https://news.ycombinator.com/item?id=9255142
#         # http://www.pixelbeat.org/programming/sigpipe_handling.html
#         set +o pipefail
#         cat {input} | sort -n -k 3 | cut -f 1 | head -n {params.top_n} \
#             | xargs kallisto index -i {output}
#         """


# rule compute_species_abundances:
#     input:
#         fastq = 'vector_depleted/{sample}.fastq',
#         index = 'kallisto_idxs/{sample}.idx'
#     output:
#         'species_abund/{sample}'
#     params:
#         length = config['insert_length'],
#         stddev = config['insert_stddev']
#     threads: 4
#     shell:
#         """
#         kallisto quant -i {input.index} -o {output} --single -l {params.length} \
#             -s {params.stddev} -t {threads} {input.fastq}
#         """
