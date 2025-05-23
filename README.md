<!-- vim-markdown-toc GFM -->

* [Description](#description)
    * [Notes on output files](#notes-on-output-files)
* [Setup & Run](#setup--run)
* [*Schistosoma* synteny tracks](#schistosoma-synteny-tracks)
* [*Trichuris* synteny tracks](#trichuris-synteny-tracks)
* [Loading assemblies, synteny, and evidence to staging server](#loading-assemblies-synteny-and-evidence-to-staging-server)
* [Local installation of Jbrowse with demo data](#local-installation-of-jbrowse-with-demo-data)
    * [Install jbrowse](#install-jbrowse)
    * [Download and add demo data](#download-and-add-demo-data)

<!-- vim-markdown-toc -->

# Description

The code here prepares demo data for showcasing
[Apollo](https://github.com/GMOD/Apollo3). You can download the processed data
from [releases](https://github.com/glaParaBio/apolloDemoData/releases).

At present there is part of the samples from [Subudhi AK et al.,
2020](https://pubmed.ncbi.nlm.nih.gov/32488076/) for *P falciparum* and *P.
chabaudi*.

## Notes on output files

(These notes may be outpdated check the Snakefile for the exact data processing)

* In directory `<ref>/hisat2`: The `bam` files include *all* reads while the
  `cram` files are downsampled to 2M reads. The zip file in releases contains
  only the cram files.

* The synteny file `crunch/blast.out.gz` is the output of `tblastx` from
  splitting the query genome in 20 kb chunks. `crunch/blast.paf` is a reshaped
  and filtered version of `blast.out.gz`. The zip file in release containbs
  only the paf file.

# Setup & Run

Install dependencies:

```
mamba create --yes -n apolloDemoData
mamba install -n apolloDemoData --yes --file requirements.txt
```

Run the analysis. Remove `-n/--dry-run` option for actual execution and remove
the `cluster` options for local execution.

```
mamba activate apolloDemoData

snakemake -p --dry-run -j 100 -C ss=$PWD/sample_sheet.tsv genomes=$PWD/genomes.tsv \
    --default-resources "mem='1G'" "cpus_per_task='1'" \
    --latency-wait 60 \
    --cluster 'sbatch -A project0014 --cpus-per-task={resources.cpus_per_task} --mem={resources.mem} --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err' \
    --cluster-cancel scancel \
    --directory ~/sharedscratch/projects/apolloDemoData \
    --keep-going
```

Once done, upload `apolloDemoData.zip` as a new release, edit tag and notes as
appropriate.

```
git add ./
git commit ...
git push
gh release create v0.4.0 ~/sharedscratch/projects/apolloDemoData/apolloDemoData.zip --notes 'Apollo demo data'
```

# *Schistosoma* synteny tracks

```
snakemake -p --dry-run -j 100 \
    -s tblastx.smk \
    --default-resources "mem='1G'" "cpus_per_task='1'" \
    --latency-wait 60 \
    --cluster 'sbatch -A project0014 --cpus-per-task={resources.cpus_per_task} --mem={resources.mem} --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err' \
    --cluster-cancel scancel \
    --directory ~/sharedscratch/projects/apolloDemoDataSchistosoma \
    --config query=schistosoma_haematobium.TD2_PRJEB44434.WBPS19 \
             subject=schistosoma_mansoni.PRJEA36577.WBPS19 \
             chroms=['1', '2', '3', '4', '5', '6', '7', 'Z', 'MITO'] \
    --keep-going
```

# *Trichuris* synteny tracks

> [!NOTE]
> Use `genomic_masked` for blast and `soft_masked` for Apollo.

```
for genome in trichuris_suis.PRJNA179528.WBPS19 trichuris_muris.PRJEB126.WBPS19 trichuris_trichiura.PRJEB535.WBPS19
do
    species=`echo ${genome} | sed 's/\..*//'`
    prj=`echo ${genome} | sed 's/.*PRJ/PRJ/ ; s/\..*//'`
    mkdir -p ${genome}/ref
    #curl https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/${species}/${prj}/${genome}.genomic_masked.fa.gz \
    #| gunzip > ${genome}/ref/${genome}.genomic_masked.fa
    #samtools faidx ${genome}/ref/${genome}.genomic_masked.fa

    curl https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/${species}/${prj}/${genome}.genomic_softmasked.fa.gz \
    | gunzip \
    | bgzip > ${genome}/ref/${genome}.genomic_softmasked.fa.gz
    samtools faidx ${genome}/ref/${genome}.genomic_softmasked.fa.gz

    ## Also prepare gff3 files. Make them smaller and remove non-CDS records with the same ID

    curl https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/${species}/${prj}/${genome}.annotations.gff3.gz \
    | zcat \
    | grep -P '\tCDS\t|\texon\t|\tgene\t|\tmRNA\t|\tprotein_coding_primary_transcript\t|\trRNA\t|\tsnRNA\t|\ttRNA\t|\tnontranslating_CDS\t|^##' > ${genome}/ref/${genome}.annotations.gff3
done
```

```
snakemake -p -j 20 \
    -s tblastx.smk \
    --default-resources "mem='1G'" "cpus_per_task='1'" \
    --latency-wait 60 \
    --cluster 'sbatch -A project0014 --cpus-per-task={resources.cpus_per_task} --mem={resources.mem} --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err' \
    --cluster-cancel scancel \
    --directory ~/sharedscratch/projects/apolloDemoDataSchistosoma \
    --config query=trichuris_muris.PRJEB126.WBPS19 \
             subject=trichuris_trichiura.PRJEB535.WBPS19 \
             chroms=['TMUE_LG2'] \
    --keep-going
```

```
snakemake -p -j 20 \
    -s tblastx.smk \
    --default-resources "mem='1G'" "cpus_per_task='1'" \
    --latency-wait 60 \
    --cluster 'sbatch -A project0014 --cpus-per-task={resources.cpus_per_task} --mem={resources.mem} --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err' \
    --cluster-cancel scancel \
    --directory ~/sharedscratch/projects/apolloDemoDataSchistosoma \
    --config query=trichuris_trichiura.PRJEB535.WBPS19 \
             subject=trichuris_suis.PRJNA179528.WBPS19 \
             chroms=['TTRE_chr2'] \
    --keep-going
```

```
snakemake -p -j 80 \
    -s tblastx.smk \
    --default-resources "mem='1G'" "cpus_per_task='1'" \
    --latency-wait 60 \
    --cluster 'sbatch -A project0014 --cpus-per-task={resources.cpus_per_task} --mem={resources.mem} --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err' \
    --cluster-cancel scancel \
    --directory ~/sharedscratch/projects/apolloDemoDataSchistosoma \
    --config query=trichuris_muris.PRJEB126.WBPS19 \
             subject=trichuris_suis.PRJNA179528.WBPS19 \
             chroms=['TMUE_LG2'] \
    --keep-going
```

# Loading assemblies, synteny, and evidence to staging server

```
cd ~/dario/manuscript
apollo login --profile demo

for fa in trichuris_muris.PRJEB126.WBPS19/ref/trichuris_muris.PRJEB126.WBPS19.genomic_softmasked.fa.gz \
    trichuris_suis.PRJNA179528.WBPS19/ref/trichuris_suis.PRJNA179528.WBPS19.genomic_softmasked.fa.gz \
    trichuris_trichiura.PRJEB535.WBPS19/ref/trichuris_trichiura.PRJEB535.WBPS19.genomic_softmasked.fa.gz
do
    name=${fa/.*/}
    apollo assembly add-from-fasta --profile demo --force ${fa} --fai ${fa}.fai --gzi ${fa}.gzi --assembly ${name}
done

apollo jbrowse get-config --profile demo > config.json

for paf in trichuris_muris.PRJEB126.WBPS19_vs_trichuris_suis.PRJNA179528.WBPS19.paf \
    trichuris_muris.PRJEB126.WBPS19_vs_trichuris_trichiura.PRJEB535.WBPS19.paf \
    trichuris_trichiura.PRJEB535.WBPS19_vs_trichuris_suis.PRJNA179528.WBPS19.paf
do
    cp ${paf} /home/ec2-user/deployment/data/
    paf=`basename $paf`
    name1=`echo $paf | sed 's/\..*//'`
    name2=`echo $paf | sed 's/.*_vs_//; s/\..*//'`

    ID1=$(
        apollo assembly get --profile demo |
        jq --raw-output ".[] | select(.name==\"$name1\")._id"
    )
    ID2=$(
        apollo assembly get --profile demo |
        jq --raw-output ".[] | select(.name==\"$name2\")._id"
    )

    jbrowse add-track \
        /data/${paf} \
      --load inPlace \
      --name "${name1} vs ${name2} TBLASTX" \
      --protocol uri \
      --assemblyNames "${ID1}","${ID2}" \
      --out config.json \
      --force
done

ID=$(
    apollo assembly get --profile demo |
    jq --raw-output ".[] | select(.name==\"trichuris_trichiura\")._id"
    )

cp trichuris_trichiura.PRJEB535.WBPS19/TTRE_all_isoseq.bam* /home/ec2-user/deployment/data/

jbrowse add-track \
  /data/TTRE_all_isoseq.bam \
  --load inPlace \
  --name "IsoSeq" \
  --protocol uri \
  --assemblyNames "${ID}" \
  --out config.json \
  --force

apollo jbrowse set-config --profile demo config.json
rm config.json
```

```
awk '$1 == "TMUE_LG2"' trichuris_muris.PRJEB126.WBPS19/ref/trichuris_muris.PRJEB126.WBPS19.annotations.gff3 > trichuris_muris.PRJEB126.WBPS19/ref/trichuris_muris.LG2.gff3
awk '$1 == "TTRE_chr2"' trichuris_trichiura.PRJEB535.WBPS19/ref/trichuris_trichiura.PRJEB535.WBPS19.annotations.gff3 > trichuris_trichiura.PRJEB535.WBPS19/ref/trichuris_trichiura.chr2.gff3

apollo feature import --profile demo trichuris_suis.PRJNA179528.WBPS19/ref/trichuris_suis.PRJNA179528.WBPS19.annotations.gff3 -a trichuris_suis -d
apollo feature import --profile demo trichuris_muris.PRJEB126.WBPS19/ref/trichuris_muris.LG2.gff3 -a trichuris_muris -d &
apollo feature import --profile demo  trichuris_trichiura.PRJEB535.WBPS19/ref/trichuris_trichiura.chr2.gff3 -a trichuris_trichiura -d &
```

# Local installation of Jbrowse with demo data

## Install jbrowse

Mostly following the jbrowse2 [quick
start](https://jbrowse.org/jb2/docs/quickstart_web/), using node v18 on Ubuntu 22 at
time of this writing:

```
npm install -g @jbrowse/cli
jbrowse create jbrowse2
cd jbrowse2/
npx serve .
```

## Download and add demo data

Inside the `jbrowse2` directory created above:

```
curl -O -L https://github.com/glaParaBio/apolloDemoData/releases/latest/download/apolloDemoData.zip
unzip -o apolloDemoData.zip
rm apolloDemoData.zip

jbrowse add-assembly apolloDemoData/Pfalciparum3D7/ref/Pfalciparum3D7.fasta --force --load copy --out demoData
jbrowse add-track apolloDemoData/Pfalciparum3D7/ref/Pfalciparum3D7.gff --force --load copy --out demoData -a Pfalciparum3D7
jbrowse add-track apolloDemoData/Pfalciparum3D7/hisat2/T0_rep1.cram --force --load copy --out demoData -a Pfalciparum3D7 --trackId T0_rep1.cram
jbrowse add-track apolloDemoData/Pfalciparum3D7/bigwig/T0_rep1.bw --force --load copy --out demoData -a Pfalciparum3D7 --trackId T0_rep1.bw

jbrowse add-assembly apolloDemoData/Pchabaudichabaudi/ref/Pchabaudichabaudi.fasta --force --load copy --out demoData
jbrowse add-track apolloDemoData/Pchabaudichabaudi/ref/Pchabaudichabaudi.gff --force --load copy --out demoData -a Pchabaudichabaudi
jbrowse add-track apolloDemoData/Pchabaudichabaudi/hisat2/Matched_D1T0830_Rep2.cram --force --load copy --out demoData -a Pchabaudichabaudi --trackId Matched_D1T0830_Rep2.cram
jbrowse add-track apolloDemoData/Pchabaudichabaudi/bigwig/Matched_D1T0830_Rep2.bw --force --load copy --out demoData -a Pchabaudichabaudi --trackId Matched_D1T0830_Rep2.bw

## Synteny file
jbrowse add-track apolloDemoData/tblastx/Pchabaudichabaudi_vs_Pfalciparum3D7.paf --assemblyNames Pchabaudichabaudi,Pfalciparum3D7 --force --load copy --out demoData

jbrowse text-index --force --out demoData
```

Open JBrowse

```
google-chrome http://localhost:3000/?config=demoData%2Fconfig.json
```
