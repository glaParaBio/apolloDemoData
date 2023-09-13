<!-- vim-markdown-toc GFM -->

* [Description](#description)
* [Setup & Run](#setup--run)

<!-- vim-markdown-toc -->

# Description

The code here prepares demo data for showcasing [Apollo](https://github.com/GMOD/Apollo3).

At present there is part of the samples from [Subudhi AK et al.,
2020](https://pubmed.ncbi.nlm.nih.gov/32488076/) for *P falciparum* and *P.
chabaudi*.

# Setup & Run

Install dependencies:

```
mamba install -n apolloDemoData --yes --file requirements.txt
```

Run the analysis. Remove `-n` (dry-run) for actual execution and remove the `cluster` options for
local execution.

```
mkdir -p output/slurm
snakemake -p -n -j 10 -C ss=$PWD/sample_sheet.tsv genomes=$PWD/genomes.tsv \
    --latency-wait 60 \
    --cluster 'sbatch --cpus-per-task=10 --mem=10G --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err' \
    --cluster-cancel scancel \
    -d output
```

Upload to github as a new release the zip directory from:

```
zip -r apolloDemoData-v0.1.zip output/demoData
```

Content:

```
├── Pchabaudichabaudi
│   ├── bigwig
│   │   └── Matched_D1T0830_Rep2.bw
│   ├── hisat2
│   │   ├── Matched_D1T0830_Rep2.cram
│   │   └── Matched_D1T0830_Rep2.cram.crai
│   └── ref
│       ├── Pchabaudichabaudi.fasta.gz
│       └── Pchabaudichabaudi.gff.gz
└── Pfalciparum3D7
    ├── bigwig
    │   └── T0_rep1.bw
    ├── hisat2
    │   ├── T0_rep1.cram
    │   └── T0_rep1.cram.crai
    └── ref
        ├── Pfalciparum3D7.fasta.gz
        └── Pfalciparum3D7.gff.gz
```
