import pandas

genomes = pandas.read_csv(config['genomes'], sep='\t', comment='#')
ss = pandas.read_csv(config['ss'], sep='\t', comment='#')
ss = ss[ss.use == 'yes']

assert len(ss.library_id) == len(set(ss.library_id))

wildcard_constraints:
    library_id='|'.join([re.escape(x) for x in ss.library_id]),
    genome='|'.join([re.escape(x) for x in ss.genome]),
    r='|'.join([re.escape(x) for x in ['1', '2']]),


rule all:
    input:
        expand('{genome}/hisat2/{library_id}.bam', zip, genome=ss.genome, library_id=ss.library_id),
        expand('{genome}/bigwig/{library_id}.bw', zip, genome=ss.genome, library_id=ss.library_id),


rule download_genome:
    output:
        fa='{genome}/ref/genome.fasta',
    params:
        url=lambda wc: genomes[genomes.genome == wc.genome].fasta.iloc[0],
    shell:
        r"""
        curl -s -L {params.url} > {output.fa}
        """


rule hisat_index:
    input:
        fa= '{genome}/ref/genome.fasta',
    output:
        idx= '{genome}/ref/genome.8.ht2',
    run:
        idx = re.sub('\.8\.ht2$', '', output.idx)

        shell(f"""
        hisat2-build -p 4 --seed 1234 -f {input.fa} {idx}
        """)



def get_fastq_download_url(ss, library_id, r):
    xs = ss[ss.library_id == library_id]
    if r == '1':
        ftp = xs.fastq_r1
    elif r == '2':
        ftp = xs.fastq_r2
    else:
        raise Exception('Invalid read mate')
    return ftp.iloc[0]


rule download_fastq:
    output:
        fq='ena/{library_id}.{srr_id}_{r}.fastq.gz',
    params:
        ftp=lambda wc: get_fastq_download_url(ss, wc.library_id, wc.r)
    shell:
        r"""
        curl -s -L {params.ftp} > {output.fq}
        """

def get_fastq_for_library_id(ss, library_id):
    fq = [ss[ss.library_id == library_id].fastq_r1.iloc[0],
          ss[ss.library_id == library_id].fastq_r2.iloc[0]]
    fq = [library_id + '.' + os.path.basename(x) for x in fq]
    fq = [os.path.join('ena', x) for x in fq]
    return fq

rule hisat2:
    input:
        fq=lambda wc: get_fastq_for_library_id(ss, wc.library_id),
        idx='{genome}/ref/genome.8.ht2',
    output:
        bam='{genome}/hisat2/{library_id}.bam',
        bai='{genome}/hisat2/{library_id}.bam.bai',
        hlog='{genome}/hisat2/{library_id}.log',
    params:
        max_intron_len= lambda wc: genomes[genomes.genome == wc.genome].max_intron_len.iloc[0] ,
    shell:
        r"""
        idx=`echo {input.idx} | sed 's/.8.ht2$//'`

        hisat2 --summary-file {output.hlog} --new-summary --fr --rna-strandness RF \
           --max-intronlen {params.max_intron_len} --threads 4 -x $idx -1 {input.fq[0]} -2 {input.fq[1]} \
        | samtools view -u -@ 4 \
        | samtools sort -@ 4 > {output.bam}
        samtools index {output.bam}
        """


rule bamToBigwig:
    input:
        bam='{genome}/hisat2/{library_id}.bam',
    output:
        bw='{genome}/bigwig/{library_id}.bw',
    shell:
        r"""
        bamCoverage -b {input.bam} -o {output} \
            --binSize 50 \
            --minMappingQuality 5 \
            --normalizeUsing BPM \
            --numberOfProcessors 4
        """
