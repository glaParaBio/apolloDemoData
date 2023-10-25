<!-- vim-markdown-toc GFM -->

* [Description](#description)
* [Setup & Run](#setup--run)
* [Local installation of Jbrowse with demo data](#local-installation-of-jbrowse-with-demo-data)
    * [Install jbrowse](#install-jbrowse)
    * [Downlaod and add demo data](#downlaod-and-add-demo-data)

<!-- vim-markdown-toc -->

# Description

The code here prepares demo data for showcasing
[Apollo](https://github.com/GMOD/Apollo3). You can download the processed data
from [releases](https://github.com/glaParaBio/apolloDemoData/releases).

At present there is part of the samples from [Subudhi AK et al.,
2020](https://pubmed.ncbi.nlm.nih.gov/32488076/) for *P falciparum* and *P.
chabaudi*.

# Setup & Run

Install dependencies:

```
mamba create --yes -n apolloDemoData
mamba install -n apolloDemoData --yes --file requirements.txt
```

Run the analysis. Remove `-n` (dry-run) for actual execution and remove the `cluster` options for
local execution.

```
mamba activate apolloDemoData

snakemake -p -n -j 10 -C ss=$PWD/sample_sheet.tsv genomes=$PWD/genomes.tsv \
    --latency-wait 60 \
    --cluster 'sbatch --cpus-per-task=10 --mem=10G --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err' \
    --cluster-cancel scancel \
    -d output
```

# Local installation of Jbrowse with demo data

## Install jbrowse

Mostly following the jbrowse2 [quick
start](https://jbrowse.org/jb2/docs/quickstart_web/), using node v18.17.1 on Ubuntu 20 at
time of this writing:

```
npm install -g @jbrowse/cli
jbrowse create jbrowse2
cd jbrowse2/
npx serve .
```

## Downlaod and add demo data

Inside the `jbrowse2` directory created above:

```
curl -O -L https://github.com/glaParaBio/apolloDemoData/releases/latest/download/apolloDemoData.zip
unzip -o apolloDemoData.zip
mv output/demoData ./
rm -r output apolloDemoData.zip

jbrowse add-assembly demoData/Pfalciparum3D7/ref/Pfalciparum3D7.fasta --load copy --out demoData
jbrowse add-track demoData/Pfalciparum3D7/ref/Pfalciparum3D7.gff.gz --load copy --out demoData -a Pfalciparum3D7
jbrowse add-track demoData/Pfalciparum3D7/hisat2/T0_rep1.cram --load copy --out demoData -a Pfalciparum3D7 --trackId T0_rep1.cram
jbrowse add-track demoData/Pfalciparum3D7/bigwig/T0_rep1.bw --load copy --out demoData -a Pfalciparum3D7 --trackId T0_rep1.bw

jbrowse add-assembly demoData/Pchabaudichabaudi/ref/Pchabaudichabaudi.fasta --load copy --out demoData
jbrowse add-track demoData/Pchabaudichabaudi/ref/Pchabaudichabaudi.gff.gz --load copy --out demoData -a Pchabaudichabaudi
jbrowse add-track demoData/Pchabaudichabaudi/hisat2/Matched_D1T0830_Rep2.cram --load copy --out demoData -a Pchabaudichabaudi --trackId Matched_D1T0830_Rep2.cram
jbrowse add-track demoData/Pchabaudichabaudi/bigwig/Matched_D1T0830_Rep2.bw --load copy --out demoData -a Pchabaudichabaudi --trackId Matched_D1T0830_Rep2.bw
```

Open JBrowse

```
google-chrome http://localhost:3000/?config=demoData%2Fconfig.json
```
