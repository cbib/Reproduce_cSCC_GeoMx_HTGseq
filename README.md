# cSCC continuum analyses

This repository contains the code required to reproduce the end-to-end analyses of **HTG EdgeSeq (bulk RNA profiling)** and **GeoMx DSP (spatial transcriptomics)** data from the publication:

> Oterino-Sogo, S. & Naji, F. et al.
> *Spatial and bulk transcriptomic profiling defines the molecular
> evolution of cutaneous squamous cell carcinoma and reveals
> stage-specific biomarkers of clinical relevance. *
<br>


# Data Availability

The analyses rely on **processed expression matrices** derived from the following GEO datasets:

-   **HTG EdgeSeq (HTG-seq)**: [GSE319969](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319969)
-   **GeoMx DSP**: [GSE319968](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319968)

⚠️ The scripts in this repository are designed to work with **processed data files**, not raw FASTQ/SRA files.
Raw sequencing files are available through the SRA links provided within each GEO entry.

<br>

# Repository Structure

The repository contains two main directories:

* `scripts/`
    * Differential expression, clustering, and visualization analyses
* `gprofiler2_ORSUM/`
    * Pathway enrichment and ORSUM summarization workflows

<br>


# Reproducing Differential Expression and Clustering Analyses

This section describes how to reproduce:

-   Differential expression analyses (LRT and Wald tests by stage; DESeq2-based)
-   LRT clustering analysis using DEGreport

## Required Input Files

1.  We provide a link to download the input files from our local storage to agilize reproducibility. You may also download the processed datasets from GEO (see **Data Availability**).
2.  Place all required input files in the directory specified by `path_to_input`.
3.  Download the file: `Heatmaps_gene_list.xlsx` and place it in the same input directory.

⚠️ The scripts will not run without these files.

## 1. Clone the Repository
```
git clone https://github.com/cbib/cSCC_continuum_analyses
cd cSCC_continuum_analyses
``` 

## 2. Download processed files

This folder contains the processed count matrices and `Heatmaps_gene_list.xlsx` file.

```
wget --no-check-certificate -r -np -nH --cut-dirs=1 -R "index.html*" http://services.cbib.u-bordeaux.fr/cSCC_gene_tables/data/
```

## 3. Create the Conda Environment
```
conda create --file geomx.yml
conda activate geomx_env
```

## 4. Configure Input and Output Paths

Edit the following scripts:

-   `GeoMx_lrt_reproduce.R`
-   `HTGseq_lrt_reproduce.R`

Update the variables:

* path_to_input
* path_to_output


## 5. Run the Analyses
```
Rscript GeoMx_lrt_reproduce.R
Rscript HTGseq_lrt_reproduce.R
```

<br>

# Output Structure

All results will be written to the directory specified in
`path_to_output` with the following structure:

* `DE_genes_by_stage/`
    * DESeq2 pairwise comparisons for HTG-seq, GeoMx Macrophages, GeoMx PanCK.
* `Heatmap_plots/`
    * Heatmaps summarizing pathway enrichment results obtained using manual curation of gprofiler2 and ORSUM output.
* `Single_gene_plots/`
    * Strip plots for candidate stage-specific cSCC biomarkers identified in the study.

Additionally, the following files are saved directly in `path_to_output`:

-   `.png` Likelihood Ratio Test (LRT) plots
-   `.csv` files containing gene groups

<br>

# Reproducing gprofiler2 and ORSUM Pathway Analysis

If you wish to reproduce gprofiler2 pathway enrichment and ORSUM summarization, follow the steps below.
⚠️ These scripts are configured to run on a **SLURM-based HPC cluster**.


## 1. Create the Environment
```
conda create --file gprofiler2_orsum.yml
conda activate gprof_orsum_env
```
## 2. Input Requirements

The script `gprofiler2_ORSUM/gprof_orsum.sh` expects a `.tsv` file as input and calls `gprofiler2_ORSUM/gprofiler_standalone.R`

Notes:

-   LRT and DEGreport outputs are saved as `.csv`. They need to be converted to `.tsv`
-   For GeoMx PanCK results, only a subset of gene groups was used for pathway analysis.
-   ORSUM requires `.gmt` files, which can be downloaded [here](https://biit.cs.ut.ee/gprofiler/gost).
-   Before running `gprof_orsum.sh` rename file prefixes by replacing `:` with `_` (example: `GO:BP` → `GO_BP`)

## 3. Run on SLURM
```
cd gprofiler2_ORSUM/

sbatch gprof_orsum.sh -i Publication_HTG_Edgeseq_gene_groups.tsv -o /results -gmt /gprofiler/gmt/hsapiens -os '500' -org 'hsapiens'

sbatch gprof_orsum.sh -i Publication_PanCK_gene_groups.tsv -o /results -gmt /gprofiler/gmt/hsapiens -os '500' -org 'hsapiens'

sbatch gprof_orsum.sh -i Publication_Macrophages_gene_groups.tsv -o /results -gmt /gprofiler/gmt/hsapiens -os '500' -org 'hsapiens'
```
Parameters:

-   `-i` : input gene group file
-   `-o` : output directory
-   `-gmt` : path to GMT files
-   `-os` : ORSUM size threshold
-   `-org` : organism (e.g., `hsapiens`)

<br>

# Contact and Support

The code in this repository was developed by
**Sergio Oterino-Sogo**

LinkedIn:
https://www.linkedin.com/in/sergio-oterino-sogo-phd-181962164/

For reproducibility issues, please open a GitHub issue:
https://github.com/cbib/cSCC_continuum_analyses/issues
