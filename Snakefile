import pandas
import glob

os.makedirs('slurm', exist_ok=True)

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
        expand('demoData/{genome}/hisat2/{library_id}.cram', zip, genome=ss.genome, library_id=ss.library_id),
        expand('demoData/{genome}/hisat2/{library_id}.cram.crai', zip, genome=ss.genome, library_id=ss.library_id),
        expand('demoData/{genome}/bigwig/{library_id}.bw', zip, genome=ss.genome, library_id=ss.library_id),
        expand('demoData/{genome}/ref/{genome}.fasta.gz', genome=ss.genome),
        expand('demoData/{genome}/ref/{genome}.gff.gz', genome=ss.genome),


rule download_genome:
    output:
        fa='{genome}/ref/{genome}.fasta',
        fai='{genome}/ref/{genome}.fasta.fai',
        dd='demoData/{genome}/ref/{genome}.fasta.gz',
    params:
        url=lambda wc: genomes[genomes.genome == wc.genome].fasta.iloc[0],
    shell:
        r"""
        curl -s -L {params.url} > {output.fa}
        gzip -c {output.fa} > {output.dd}
        samtools faidx {output.fa}
        """


rule download_gff:
    output:
        gff='demoData/{genome}/ref/{genome}.gff.gz',
    params:
        url=lambda wc: genomes[genomes.genome == wc.genome].gff.iloc[0],
    shell:
        r"""
        curl -s -L {params.url} | gzip > {output.gff}
        """


rule hisat_index:
    input:
        fa= '{genome}/ref/{genome}.fasta',
    output:
        idx= '{genome}/ref/{genome}.8.ht2',
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
        idx='{genome}/ref/{genome}.8.ht2',
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


rule cram_for_demo:
    input:
        bam='{genome}/hisat2/{library_id}.bam',
        genome='{genome}/ref/{genome}.fasta',
    output:
        cram='demoData/{genome}/hisat2/{library_id}.cram',
        crai='demoData/{genome}/hisat2/{library_id}.cram.crai',
    params:
        ctg=lambda wc: genomes[genomes.genome== wc.genome].demo_contigs.iloc[0].replace(',', ' '),
    shell:
        r"""
        samtools view -C -T {input.genome} {input.bam} {params.ctg} > {output.cram}
        samtools index {output.cram}
        """


rule bamToBigwig:
    input:
        bam='{genome}/hisat2/{library_id}.bam',
    output:
        bw='demoData/{genome}/bigwig/{library_id}.bw',
    shell:
        r"""
        bamCoverage -b {input.bam} -o {output} \
            --binSize 50 \
            --minMappingQuality 5 \
            --normalizeUsing BPM \
            --numberOfProcessors 4
        """


rule makeblastdb:
    input:
        fa='Pfalciparum3D7/ref/Pfalciparum3D7.fasta',
    output:
        fa='blast/db/ref.fasta',
        db='blast/db/ref.fasta.nin',
    params:
        ctg=lambda wc: genomes[genomes.genome== 'Pfalciparum3D7'].demo_contigs.iloc[0].replace(',', ' '),
    shell:
        r"""
        samtools faidx {input.fa} {params.ctg} > {output.fa}
        makeblastdb -in {output.fa} -dbtype nucl
        """


checkpoint split_query:
    input:
        fai='Pchabaudichabaudi/ref/Pchabaudichabaudi.fasta.fai',
        fa='Pchabaudichabaudi/ref/Pchabaudichabaudi.fasta',
    output:
        outdir=temp(directory('split')),
    params:
        ctg=lambda wc: genomes[genomes.genome== 'Pchabaudichabaudi'].demo_contigs.iloc[0].replace(',', '|'),
    shell:
        r"""
        mkdir {output.outdir}
        grep -w -P '{params.ctg}' {input.fai} | bedtools makewindows -g - -w 100000 | awk '{{print $1 ":" $2+1 "-" $3}}' > {output.outdir}/windows.bed
        while read -r line
        do
           samtools faidx {input.fa} $line > split/${{line}}.fasta
        done < {output.outdir}/windows.bed
        rm {output.outdir}/windows.bed
        """


rule tblastx:
    input:
        query='split/{window}.fasta',
        db='blast/db/ref.fasta.nin',
        fa='blast/db/ref.fasta',
    output:
        out='blast/{window}.out',
    shell:
        r"""
        tblastx -query {input.query} \
           -db {input.fa} \
           -evalue 0.1 \
           -max_target_seqs 10 \
           -outfmt 6 > {output.out} 
        rm {input.query}
        """


def aggregate_blast(wc):
    checkpoint_output = checkpoints.split_query.get().output.outdir
    outfiles = [os.path.basename(x) for x in glob.glob(os.path.join(checkpoint_output, '*.fasta'))]
    windows = [re.sub('\.fasta$', '', x) for x in outfiles]
    return expand('blast/{window}.out', window=windows)

rule cat_blast:
    input:
        out=aggregate_blast,
    output:
        out='crunch/blast.out.gz',
    run:
        fout = open(output.out + '.tmp', 'w')
        for xin in input.out:
            with open(xin) as fin:
                for line in fin:
                    if line.startswith('#'):
                        continue
                    line = line.strip().split('\t')
                    ctg = line[0].split(':')
                    assert len(ctg) == 2
                    ctg = ctg[0]
                    offset = re.sub('-.*', '', re.sub('.*:', '', line[0]))
                    offset = int(offset) - 1
                    assert offset >= 0
                    line[0] = ctg
                    line[6] = str(int(line[6])  + offset)
                    line[7] = str(int(line[7])  + offset)
                    fout.write('\t'.join(line) + '\n')
        fout.close()
        shell('sort -k1,1 -k7,7n -k9,9n {output.out}.tmp | gzip > {output.out} && rm {output.out}.tmp')
