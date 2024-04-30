<!-- vim-markdown-toc GFM -->

* [Description](#description)
    * [Notes on output files](#notes-on-output-files)
* [Setup & Run](#setup--run)
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

snakemake -p --dry-run -j 500 -C ss=$PWD/sample_sheet.tsv genomes=$PWD/genomes.tsv \
    --default-resources "mem='1G'" "cpus_per_task='4'" \
    --latency-wait 60 \
    --cluster 'sbatch --cpus-per-task={resources.cpus_per_task} --mem={resources.mem} --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err' \
    --cluster-cancel scancel \
    --directory ~/sharedscratch/projects/apolloDemoData
```

Once done, upload `apolloDemoData.zip` as a new release, edit tag and notes as
appropriate.

```
git add ./
git commit ...
git push
gh release create v0.4.0 ~/sharedscratch/projects/apolloDemoData/apolloDemoData.zip --notes 'Apollo demo data'
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
