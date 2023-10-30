import pandas
import glob

os.makedirs('slurm', exist_ok=True)

genomes = pandas.read_csv(config['genomes'], sep='\t', comment='#')
ss = pandas.read_csv(config['ss'], sep='\t', comment='#')
ss = ss[ss.use == 'yes']

BLAST_FMT = 'qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore nident positive slen sstrand'

assert len(ss.library_id) == len(set(ss.library_id))

wildcard_constraints:
    library_id='|'.join([re.escape(x) for x in ss.library_id]),
    genome='|'.join([re.escape(x) for x in ss.genome]),
    r='|'.join([re.escape(x) for x in ['1', '2']]),

localrules: all, download_fastq, download_gff

rule all:
    input:
        'apolloDemoData.zip',
        'apolloDemoDataFull.zip',


rule zip_demodata:
    input:
        expand('{genome}/hisat2/{library_id}.cram', zip, genome=ss.genome, library_id=ss.library_id),
        expand('{genome}/hisat2/{library_id}.cram.crai', zip, genome=ss.genome, library_id=ss.library_id),
        expand('{genome}/bigwig/{library_id}.bw', zip, genome=ss.genome, library_id=ss.library_id),
        expand('{genome}/ref/{genome}.fasta', genome=ss.genome),
        expand('{genome}/ref/{genome}.fasta.fai', genome=ss.genome),
        expand('{genome}/ref/{genome}.gff.gz', genome=ss.genome),
        expand('{genome}/ref/{genome}.gff.gz.tbi', genome=ss.genome),
        'crunch/blast.paf',
    output:
        zip='apolloDemoData.zip',
    params:
        outdir=lambda wc, output: re.sub('\.zip$', '', output.zip)
    shell:
        r"""
        rm -rf {params.outdir}
        mkdir {params.outdir}
        rsync --relative -arvP {input} {params.outdir}
        zip -r {output.zip} {params.outdir}
        rm -r {params.outdir}
        """


rule zip_fulldemodata:
    input:
        expand('{genome}/hisat2/{library_id}.cram', zip, genome=ss.genome, library_id=ss.library_id),
        expand('{genome}/hisat2/{library_id}.cram.crai', zip, genome=ss.genome, library_id=ss.library_id),
        expand('{genome}/hisat2/{library_id}.bam', zip, genome=ss.genome, library_id=ss.library_id),
        expand('{genome}/hisat2/{library_id}.bam.bai', zip, genome=ss.genome, library_id=ss.library_id),
        expand('{genome}/bigwig/{library_id}.bw', zip, genome=ss.genome, library_id=ss.library_id),
        expand('{genome}/ref/{genome}.fasta', genome=ss.genome),
        expand('{genome}/ref/{genome}.fasta.fai', genome=ss.genome),
        expand('{genome}/ref/{genome}.gff.gz', genome=ss.genome),
        expand('{genome}/ref/{genome}.gff.gz.tbi', genome=ss.genome),
        'crunch/blast.paf',
        'crunch/blast.out.gz',
    output:
        zip='apolloDemoDataFull.zip',
    params:
        outdir=lambda wc, output: re.sub('\.zip$', '', output.zip)
    resources:
        mem='2G',
    shell:
        r"""
        rm -rf {params.outdir}
        mkdir {params.outdir}
        rsync --relative -arvP {input} {params.outdir}
        zip -r {output.zip} {params.outdir}
        rm -r {params.outdir}
        """


rule download_genome:
    output:
        fa='{genome}/ref/{genome}.fasta',
        fai='{genome}/ref/{genome}.fasta.fai',
    params:
        url=lambda wc: genomes[genomes.genome == wc.genome].fasta.iloc[0],
    resources:
        mem='100M',
    shell:
        r"""
        curl -s -L {params.url} > {output.fa}
        samtools faidx {output.fa}
        """


rule download_gff:
    output:
        gff='{genome}/ref/{genome}.gff.gz',
    params:
        url=lambda wc: genomes[genomes.genome == wc.genome].gff.iloc[0],
    resources:
        mem='100M',
    shell:
        r"""
        curl -s -L {params.url} \
        | sortBed -header \
        | bgzip > {output.gff}
        """


rule index_gff:
    input:
        gff='{genome}/ref/{genome}.gff.gz',
    output:
        tbi='{genome}/ref/{genome}.gff.gz.tbi',
    shell:
        r"""
        tabix -p gff {input.gff}
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
    resources:
        mem='100M',
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
    resources:
        mem='6G',
    shell:
        r"""
        idx=`echo {input.idx} | sed 's/.8.ht2$//'`

        hisat2 --summary-file {output.hlog} --new-summary --fr --rna-strandness RF \
           --max-intronlen {params.max_intron_len} --threads 4 -x $idx -1 {input.fq[0]} -2 {input.fq[1]} \
        | samtools view -u -@ 4 \
        | samtools sort -@ 4 > {output.bam}
        samtools index {output.bam}
        """


rule downsample_cram:
    input:
        bam='{genome}/hisat2/{library_id}.bam',
        genome='{genome}/ref/{genome}.fasta',
    output:
        cram='{genome}/hisat2/{library_id}.cram',
    shell:
        r"""
        set +o pipefail
        tot=`samtools stats -@ 4 {input.bam} | grep '1st fragments:' | cut -f 3`
        set -o pipefail
        frac=`echo 2000000 / $tot | bc -l | awk '{{if($1 > 1){{print 1}} else {{print $1}}}}'`
         
        samtools view -@ 4 --subsample $frac --subsample-seed 1234 -C -T {input.genome} {input.bam} > {output.cram}
        """


rule index_cram:
    input:
        cram='{genome}/hisat2/{library_id}.cram',
    output:
        crai='{genome}/hisat2/{library_id}.cram.crai',
    shell:
        r"""
        samtools index {input.cram}
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


rule makeblastdb:
    input:
        fa='Pfalciparum3D7/ref/Pfalciparum3D7.fasta',
    output:
        db='Pfalciparum3D7/ref/Pfalciparum3D7.fasta.nin',
    #params:
    #    ctg=lambda wc: genomes[genomes.genome== 'Pfalciparum3D7'].demo_contigs.iloc[0].replace(',', ' '),
    shell:
        r"""
        makeblastdb -in {input.fa} -dbtype nucl
        """


checkpoint split_query:
    input:
        fai='Pchabaudichabaudi/ref/Pchabaudichabaudi.fasta.fai',
        fa='Pchabaudichabaudi/ref/Pchabaudichabaudi.fasta',
    output:
        outdir=temp(directory('split')),
    #params:
    #    ctg=lambda wc: genomes[genomes.genome== 'Pchabaudichabaudi'].demo_contigs.iloc[0].replace(',', '|'),
    shell:
        r"""
        mkdir {output.outdir}
        bedtools makewindows -g {input.fai} -w 20000 | awk '{{print $1 ":" $2+1 "-" $3}}' > {output.outdir}/windows.bed
        while read -r line
        do
           samtools faidx {input.fa} $line > split/${{line}}.fasta
        done < {output.outdir}/windows.bed
        rm {output.outdir}/windows.bed
        """


rule tblastx:
    input:
        query='split/{window}.fasta',
        db='Pfalciparum3D7/ref/Pfalciparum3D7.fasta.nin',
        fa='Pfalciparum3D7/ref/Pfalciparum3D7.fasta',
    output:
        out='blast/{window}.out',
    params:
        fmt=BLAST_FMT,
    resources:
        mem='1G',
        cpus_per_task=1,
    shell:
        r"""
        tblastx -query {input.query} \
           -db {input.fa} \
           -evalue 0.1 \
           -max_target_seqs 10 \
           -outfmt "6 {params.fmt}" > {output.out} 
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
    resources:
        mem='2500M',
    run:
        import gzip 
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
        
        HEADER = BLAST_FMT.replace(' ', '\t')
        with gzip.open(output.out, 'wb') as fout:
            fout.write(f'#{HEADER}\n'.encode())

        shell('LC_ALL=C sort -S 2G --parallel 4 -k1,1 -k7,7n -k9,9n {output.out}.tmp | gzip >> {output.out} && rm {output.out}.tmp')


rule blast_to_paf:
    input:
        out='crunch/blast.out.gz',
        fai=expand('{genome}/ref/{genome}.fasta.fai', genome=ss.genome),
    output:
        paf='crunch/blast.paf',
    resources:
        mem='8G',
    run:
        fai = {}
        for fn in input.fai:
            with open(fn) as fin:
                for line in fin:
                    line = line.strip().split('\t')
                    assert line[0] not in fai
                    fai[line[0]] = line[1]
        
        out = pandas.read_csv(input.out, sep='\t')
        out = out[(out.pident >= 50) & (out.evalue < 0.01) & (out.length >= 20)]
        out.rename(columns={'#qaccver': 'qaccver'}, inplace=True) 

        qstart = []
        qend = []
        sstart = []
        send = []
        qaccver_length = []
        saccver_length = []
        strand = []
        for row in out.itertuples():
            qaccver_length.append(fai[row.qaccver])
            saccver_length.append(fai[row.saccver])
            if ((row.qstart < row.qend) and (row.sstart < row.send)) or ((row.qstart > row.qend) and (row.sstart > row.send)):
                strand.append('+')
            else:
                strand.append('-')

            if row.qstart < row.qend:
                qstart.append(row.qstart)
                qend.append(row.qend)
            else:
                qstart.append(row.qend)
                qend.append(row.qstart)
            if row.sstart < row.send:
                sstart.append(row.sstart)
                send.append(row.send)
            else:
                sstart.append(row.send)
                send.append(row.sstart)
        out['qstart'] = [x - 1 for x in qstart]
        out['sstart'] = [x - 1 for x in sstart]
        out['qend'] = qend
        out['send'] = send
        out['qaccver_length'] = qaccver_length
        out['saccver_length'] = saccver_length
        out['strand'] = strand
        out['aln_length'] = out['qend'] - out['qstart']
        out['mapq'] = 255

        out.rename(columns={'qaccver': '#qaccver'}, inplace=True)
        out[['#qaccver', 'qaccver_length', 'qstart', 'qend', 'strand', 'saccver', 'saccver_length', 'sstart', 'send', 'nident', 'aln_length', 'mapq']].to_csv(output.paf, sep='\t', index=False)

