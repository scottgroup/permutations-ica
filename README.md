## Factorial study of the RNA-seq computational workflow identifies biases as technical gene signatures

This repository contains the Snakemake project of the following paper :

> Joël Simoneau, Ryan Gosselin, Michelle S. Scott.
> "Factorial study of the RNA-seq computational workflow identifies biases as technical gene signatures" *bioRxiv*, January 2020.
> [doi: 10.1101/2020.01.30.924092](https://doi.org/10.1101/2020.01.30.924092)

### How to run
This project is wrapped as a Snakemake workflow using conda to manage the different software environments needed for the analyses.

1 - Installing conda (for Linux users) :
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Answer `yes` to `Do you wish the installer to initialize Miniconda3?`


2 - Creating a Snakemake environment. To create the Snakemake environment used to launch Snakemake, run the following. The `conda create` command can appear to be stuck on `Solving environment`. While we are actually arguably [never going to solve the environment](https://www.ipcc.ch/sr15/chapter/spm/), the command is probably not stuck. Just be patient.

```bash
exec bash
conda config --set auto_activate_base False
conda create --name smake -c bioconda -c conda-forge snakemake=5.7.4
```

Before running Snakemake, you have to initialize the environment
```bash
conda activate smake
```

3 - Run the project. To run the workflow locally simply run the following command in the Snakemake conda environment, where `$CORES` is the number of available cores.
```bash
snakemake --use-conda --cores=$CORES
```

To run on a Slurm cluster, one can use the following command to output all tasks at once, where `$MAXNJOBS` represents that maximum number of jobs a user can launch at the same time.
```bash
snakemake -j $MAXNJOBS --use-conda --immediate-submit --notemp --cluster-config subworkflows/rna_seq_cartesian_product/cluster.json --cluster-config cluster.json  --cluster 'python3 slurmSubmit.py {dependencies}'
```

The `cluster.json` files need to be ajusted accordingly to the cluster architecture.

### Structure of the project
The main Snakemake workflow describes the generation and analysis of the ICA models. This workflow contains two subworkflows. The first one is responsible for the cartesian product processing of the RNA-seq pipelines (`subworkflows/rna_seq_cartesian_product/`) and the second generates the Supplementary Data 3 (`subworkflows/refseq_noexon/`).

#### RNA-seq cartesian product
This subworkflow downloads the datasets, process them and summarizes the counts in a file that is passed to the main workflow.

The following software and references were included in this study :
1. Trimming: Cuptadapt, Trimmomatic
2. Genome annotation: Ensembl 92, Ensembl 98 and RefSeq
3. Alignment: HISAT2, STAR, TopHat2
4. Quantification: Cufflinks, featureCounts, HTSeq

This subworkflow generates two different files that can be used by the ICA analysis.

- `subworkflows/rna_seq_cartesian_product/results/cartesian_product/tissues_NaN.tsv` contains the gene-level read counts for genes that were reported in all the different pipelines.

- `subworkflows/rna_seq_cartesian_product/results/cartesian_product/tissues_noNaN.tsv` contains the gene-level read counts for genes that were reported in at least one pipeline. Missing values are filled with 0.

#### ICA models
The main workflow generates the ICA models, and provides the tools to analyse them.

```javascript
"ICA_models":
{
    "tissues_NaN_tripleAnnot_cufflinks":
    {
        "params":
        {
            "counts": "tissues_NaN",
            "ICAmethod": ["sklearnFastICA"],
            "max": 35,
            "min": 6,
            "n": 25,
            "sigma": [4],
            "std": [0],
            "threshold": 0.90
        },
        "variables":
        {
            "annotation": ["ensembl92", "ensembl98", "refseq"],
            "quantifier": ["cufflinks"]
        }
    }
}
```
