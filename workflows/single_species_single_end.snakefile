from glob import glob

from llutil.utils import fastx_stem


config['input_glob'] =
config['insert_length'] =
config['species_index'] =
config['species_gff'] =
config['species_fasta'] =
config['read_length'] =


input_files = glob(config['input_glob'])


rule all:
    input:
        'genome_aligned_reads.unsorted.bam',
        'fragments.tsv',
        'coverage/genome.weighted.tsv',
        'coverage/genome.dedup.tsv',
        'coverage/unstranded.weighted.tsv',
        'coverage/unstranded.dedup.tsv',
        'coverage/stranded.weighted.tsv',
        'coverage/stranded.dedup.tsv',
        'coverage/inframe.weighted.tsv',
        'coverage/inframe.dedup.tsv',
        # 'clusters.counts.tsv',
        # 'annot/{sample}.inframe_inserts.counts.tsv', sample=config['samples']),
        # 'ref/{sample}.uniq_frag.bed', sample=config['samples']),
        # expand('diversity/{sample}.recon_allreads.txt', sample=config['samples']),
        # expand('diversity/{sample}.recon_inframe.txt', sample=config['samples'])


# GENERATE REFERENCE DATA


rule cds_filter_gff:
    input:
        config['species_gff']
    output:
        temp('cdss.gff')
    shell:
        r"""cat {input} | awk -F \t '/^[^#]/ && $3 == "CDS"' > {output}"""


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


# COMPUTE ALIGNMENTS AND ANNOTATED FRAGMENT COUNTS


rule align_reads_to_genome:
    input:
        fastq = input_files
    output:
        'genome_aligned_reads.unsorted.bam'
    params:
        species_index = config['species_index']
    threads: 4
    shell:
        """
        cat {input} \
            | bowtie2 -p {threads} -x {params.species_index} -U - \
            | samtools view -b - \
            > {output}
        """


rule generate_fragments:
    input:
        alignments = 'genome_aligned_reads.unsorted.bam',
        genome = 'chrom_sizes.tsv',
        cds_annot = 'cdss.gff'
    output:
        temp('raw_fragments.tsv')
    params:
        size = config['insert_length']
    shell:
        r"""
        echo 'UNDO $3 != 0 when bedtools2#646 is fixed'
        echo -e 'chr\tfrag_start\tfrag_end\tstrand\tcds_start\tcds_end\tcount\tannot' > {output}
        bedtools bamtobed -i {input.alignments} \
            | awk -F '\t' 'BEGIN {{OFS="\t"}} {{print $1, ($6 == "+") ? $2 : $3, $6}}' \
            | sort \
            | uniq -c \
            | awk 'BEGIN {{OFS="\t"}} $3 != 0 {{print $2, $3, $3, ".", $1, $4}}' \
            | bedtools slop -s -l 0 -r {params.size} -i stdin -g {input.genome} \
            | bedtools intersect -s -wao -a stdin -b {input.cds_annot} \
            | awk -F '\t' 'BEGIN {{OFS="\t"}} {{print $1, $2, $3, $6, $10, $11, $5, $15}}' \
            >> {output}
        """


rule annotate_fragments:
    input:
        'raw_fragments.tsv',
        config['species_fasta']
    output:
        'fragments.tsv'
    run:
        from csv import DictReader
        from collections import Counter

        from tqdm import tqdm
        from Bio import SeqIO

        chrom_dict = SeqIO.to_dict(SeqIO.parse(input[1], 'fasta'))

        fieldnames = ['chr', 'frag_start', 'frag_end', 'strand', 'cds_start', 'cds_end', 'count', 'annot']
        annots = dict()
        with open(input[0], 'r') as ip, open(output[0], 'w') as op:
            _ = next(ip)  # skip input header
            # write output header
            print('chr\tposition\tstrand\tcount\thas_overlap\thas_inframe_overlap\tbases_till_overlap\tfragment_nt\tfragment_aa\tannot', file=op)
            for r in tqdm(DictReader(ip, fieldnames=fieldnames, delimiter='\t', dialect='unix')):
                if r['strand'] not in ['+', '-']:
                    raise ValueError('strand must be + or -')

                cds_start = int(r['cds_start'])
                cds_end = int(r['cds_end'])
                frag_start = int(r['frag_start'])
                frag_end = int(r['frag_end'])

                # is the fragment on the positive strand of the chr?
                is_positive_strand = int(r['strand'] == '+')
                # does the fragment overlap any features in the bed file?
                # note: this value comes from bedtools, which uses -1 for no
                # overlap
                has_overlap = int(cds_start >= 0)
                # the python-style index where the fragment starts
                position = frag_start if is_positive_strand else frag_end

                # compute base differential between fragment and CDS
                if is_positive_strand:
                    # off-by-one bc 'position' derives from a BED file and
                    # 'cds_start' derives from a GFF file
                    diff = position - cds_start + 1
                else:
                    diff = cds_end - position

                # extract nt and aa sequences
                fragment_nt = chrom_dict[r['chr']][frag_start:frag_end].seq
                if not is_positive_strand:
                    fragment_nt = fragment_nt.reverse_complement()
                frame = len(fragment_nt) % 3
                fragment_aa = fragment_nt[:len(fragment_nt) - frame].translate().split('*')[0]

                # how many bases from the 5' end of the fragment till we
                # get to the overlap
                bases_till_overlap = max(-diff, 0) if has_overlap else ''
                correct_frame = diff % 3 == 0
                translation_reaches_overlap = (len(fragment_aa) * 3 > bases_till_overlap) if has_overlap else False
                has_inframe_overlap = int(has_overlap and correct_frame and translation_reaches_overlap)

                record = (
                    r['chr'], position, r['strand'], r['count'], has_overlap,
                    has_inframe_overlap, bases_till_overlap, fragment_nt,
                    fragment_aa, r['annot'])
                record_string = '\t'.join(map(str, record))
                print(record_string, file=op)


# COMPUTE COVERAGE


rule genome_coverage:
    input:
        fragments = 'raw_fragments.tsv',
        genome = 'chrom_sizes.tsv'
    output:
        weighted = 'coverage/genome.weighted.tsv',
        dedup = 'coverage/genome.dedup.tsv'
    shell:
        """
        tail -n +2 {input.fragments} \
            | awk -F '\t' 'BEGIN {{OFS="\t"}} {{for (i=0; i < $7; i++) print $1, $2, $3, ".", ".", $4}}' \
            | bedtools genomecov -d -g {input.genome} -i stdin \
            > {output.weighted}
        tail -n +2 {input.fragments} \
            | awk -F '\t' 'BEGIN {{OFS="\t"}} {{print $1, $2, $3, ".", ".", $4}}' \
            | bedtools genomecov -d -g {input.genome} -i stdin \
            > {output.dedup}
        """


rule unstranded_cds_coverage:
    input:
        fragments = 'raw_fragments.tsv',
        cdss = 'cdss.gff'
    output:
        weighted = 'coverage/unstranded.weighted.tsv',
        dedup = 'coverage/unstranded.dedup.tsv'
    shell:
        """
        tail -n +2 {input.fragments} \
            | awk -F '\t' 'BEGIN {{OFS="\t"}} {{for (i=0; i < $7; i++) print $1, $2, $3, ".", ".", $4}}' \
            | bedtools coverage -a {input.cdss} -b stdin -wao \
            > {output.weighted}
        tail -n +2 {input.fragments} \
            | awk -F '\t' 'BEGIN {{OFS="\t"}} {{print $1, $2, $3, ".", ".", $4}}' \
            | bedtools coverage -a {input.cdss} -b stdin -wao \
            > {output.dedup}
        """


rule stranded_cds_coverage:
    input:
        fragments = 'raw_fragments.tsv',
        cdss = 'cdss.gff'
    output:
        weighted = 'coverage/stranded.weighted.tsv',
        dedup = 'coverage/stranded.dedup.tsv'
    shell:
        """
        tail -n +2 {input.fragments} \
            | awk -F '\t' 'BEGIN {{OFS="\t"}} {{for (i=0; i < $7; i++) print $1, $2, $3, ".", ".", $4}}' \
            | bedtools coverage -s -a {input.cdss} -b stdin -wao \
            > {output.weighted}
        tail -n +2 {input.fragments} \
            | awk -F '\t' 'BEGIN {{OFS="\t"}} {{print $1, $2, $3, ".", ".", $4}}' \
            | bedtools coverage -s -a {input.cdss} -b stdin -wao \
            > {output.dedup}
        """


rule inframe_cds_coverage:
    input:
        fragments = 'raw_fragments.tsv',
        cdss = 'cdss.gff'
    output:
        weighted = 'coverage/inframe.weighted.tsv',
        dedup = 'coverage/inframe.dedup.tsv'
    shell:
        """
        tail -n +2 {input.fragments} \
            | awk '$5 >= 0' \
            | awk '($4 == "+" && ($2 - $5 + 1) % 3 == 0) || ($4 == "-" && ($6 - $3) % 3 == 0)' \
            | awk -F '\t' 'BEGIN {{OFS="\t"}} {{for (i=0; i < $7; i++) print $1, $2, $3, ".", ".", $4}}' \
            | bedtools coverage -s -a {input.cdss} -b stdin -wao \
            > {output.weighted}
        tail -n +2 {input.fragments} \
            | awk '$5 >= 0' \
            | awk '($4 == "+" && ($2 - $5 + 1) % 3 == 0) || ($4 == "-" && ($6 - $3) % 3 == 0)' \
            | awk -F '\t' 'BEGIN {{OFS="\t"}} {{print $1, $2, $3, ".", ".", $4}}' \
            | bedtools coverage -s -a {input.cdss} -b stdin -wao \
            > {output.dedup}
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


# rule recon_diversity_allreads:
#     input:
#         'clustered_reads/{sample}.clustered.counts.tsv'
#     output:
#         temp('diversity/tmp_{sample}.recon_input.tsv'),
#         'diversity/{sample}.recon_allreads.txt'
#     params:
#         python2 = config['python2'],
#         recon = config['recon']
#     shell:
#         r"""
#         cat {input} | awk '{{print $2 "\t" $1}}' > {output[0]}
#         {params.python2} {params.recon} -R -t 1000 -o {output[1]} {output[0]}
#         """


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


# rule recon_diversity_inframe:
#     input:
#         'annot/{sample}.inframe_inserts.counts.tsv'
#     output:
#         'diversity/{sample}.recon_inframe.txt'
#     params:
#         python2 = config['python2'],
#         recon = config['recon']
#     shell:
#         r"""
#         {params.python2} {params.recon} -R -t 1000 -o {output} {input}
#         """

rule gen_align_ref:
    input:
        '{sample}.annot_inserts.tsv'
    output:
        'ref/{sample}.uniq_frag.bed'
    shell:
        r"""
        # First awk script pulls out only relevant columns (and computes
        # in-frame) into bed6 format:
        # chrom, start, end, annot, inframe, strand
        # Second awk script output only uniq lines (no sort req'd)
        cat {input} \
            | awk -F '\t' 'BEGIN {{OFS="\t"}} {{inframe = ($11 >= 0) && (($6 == "+" && ($2 - $11 + 1) % 3 == 0) || ($6 == "-" && ($3 - $12) % 3 == 0)); print $1, $2, $3, $16, inframe, $6}}' \
            | awk '!($0 in seen) {{++seen[$0]; print $0}}' \
            > {output}
        """
