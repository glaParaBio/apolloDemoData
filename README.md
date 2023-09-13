
```
mamba install -n apolloDemoData --yes --file requirements.txt
```

# Run

```
mkdir -p output/slurm
snakemake -p -n -j 10 -C ss=$PWD/sample_sheet.tsv genomes=$PWD/genomes.tsv \
    --latency-wait 60 \
    --cluster 'sbatch --cpus-per-task=10 --mem=10G --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err' \
    --cluster-cancel scancel \
    -d output
```
