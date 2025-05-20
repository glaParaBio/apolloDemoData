import pandas
import glob

os.makedirs('slurm', exist_ok=True)

CHROMS=config['chroms']

print(CHROMS)

QUERY=config['query'] 
SUBJECT=config['subject'] 

BLAST_FMT = 'qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore nident positive slen sstrand'

wildcard_constraints:
    chrom='|'.join([re.escape(x) for x in CHROMS]),

localrules: all

rule all:
    input:
        expand('tblastx/{chrom}/{QUERY}_vs_{SUBJECT}.paf', QUERY=QUERY, SUBJECT=SUBJECT, chrom=CHROMS),

rule faidx:
    input:
        genome='{genome}/ref/{genome}.genomic_masked.fa',
    output:
        fai='{genome}/ref/{genome}.genomic_masked.fa.fai',
    shell:
        r"""
        samtools faidx {input}
        """

rule getChrom:
    input:
        genome='{genome}/ref/{genome}.genomic_masked.fa',
        fai='{genome}/ref/{genome}.genomic_masked.fa.fai',
    output:
        chrom=temp('{genome}/ref/{chrom}.fa'),
        fai=temp('{genome}/ref/{chrom}.fa.fai'),
    shell:
        r"""
        chrom=`grep -P '^{wildcards.chrom}\t' {input.fai} | cut -f 1`
        samtools faidx {input.genome} $chrom > {output.chrom}
        samtools faidx {output.chrom}
        """

rule makeblastdb:
    input:
        fa='{subject}/ref/{subject}.genomic_masked.fa',
        # fa='{subject}/ref/{chrom}.fa',
    output:
        db='{subject}/ref/{subject}.genomic_masked.fa.nin',
        # db='{subject}/ref/{chrom}.fa.nin',
    shell:
        r"""
        makeblastdb -in {input.fa} -dbtype nucl
        """


checkpoint split_query:
    input:
        fai='{query}/ref/{chrom}.fa.fai',
        fa='{query}/ref/{chrom}.fa',
    output:
        outdir=temp(directory('split/{query}/{chrom}')),
    shell:
        r"""
        mkdir {output.outdir}
        bedtools makewindows -g {input.fai} -w 20000 | awk '{{print $1 ":" $2+1 "-" $3}}' > {output.outdir}/windows.bed
        while read -r line
        do
           samtools faidx {input.fa} $line > {output.outdir}/${{line}}.fa
        done < {output.outdir}/windows.bed
        rm {output.outdir}/windows.bed
        """


rule tblastx:
    input:
        query='split/{query}/{chrom}/{window}.fa',
        db='{subject}/ref/{subject}.genomic_masked.fa.nin',
        fa='{subject}/ref/{subject}.genomic_masked.fa',
        # db='{subject}/ref/{chrom}.fa.nin',
        # fa='{subject}/ref/{chrom}.fa',
    output:
        out='blast/{chrom}/{query}/{subject}/{window}.out',
    params:
        fmt=BLAST_FMT,
    resources:
        mem='4G',
        cpus_per_task=1,
        time='04:00:00',
    shell:
        r"""
        tblastx -query {input.query} \
           -db {input.fa} \
           -evalue 0.1 \
           -max_target_seqs 10 \
           -outfmt "6 {params.fmt}" > {output.out} 
        """


def aggregate_blast(wc):
    checkpoint_output = checkpoints.split_query.get(query=wc.query, chrom=wc.chrom).output.outdir
    outfiles = [os.path.basename(x) for x in glob.glob(os.path.join(checkpoint_output, '*.fa'))]
    windows = [re.sub('\.fa$', '', x) for x in outfiles]
    return expand('blast/{chrom}/{query}/{subject}/{window}.out', chrom=wc.chrom, window=windows, query=wc.query, subject=wc.subject)

rule cat_blast:
    input:
        out=aggregate_blast,
    output:
        out='tblastx/{chrom}/{query}_vs_{subject}.out.gz',
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
        out='tblastx/{chrom}/{query}_vs_{subject}.out.gz',
        fai_query='{query}/ref/{chrom}.fa.fai',
        fai_subject='{subject}/ref/{subject}.genomic_masked.fa.fai',
    output:
        paf='tblastx/{chrom}/{query}_vs_{subject}.paf',
    resources:
        mem='8G',
    run:
        fai = {}
        for fn in [input.fai_query, input.fai_subject]:
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

