# cSCC Continuum Analyses on HTG EdgeSeq and GeoMx DSP Data

This repository contains the code required to reproduce the end-to-end analyses of **HTG EdgeSeq (bulk RNA profiling)** and **GeoMx DSP (spatial transcriptomics)** data from the publication:

> Oterino-Sogo, S. & Naji, F. et al.\
> *Spatial and bulk transcriptomic profiling defines the molecular
> evolution of cutaneous squamous cell carcinoma and reveals
> stage-specific biomarkers of clinical relevance.*

------------------------------------------------------------------------

# Data Availability

The analyses rely on **processed expression matrices** derived from the following GEO datasets:

-   **HTG EdgeSeq (HTG-seq)**:\
    GSE319969\
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319969

-   **GeoMx DSP**:\
    GSE319968\
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319968

⚠️ The scripts in this repository are designed to work with **processed data files**, not raw FASTQ/SRA files.\
Raw sequencing files are available through the SRA links provided within each GEO entry.

------------------------------------------------------------------------

# Repository Structure

The repository contains two main directories:

* scripts/
* gprofiler2_ORSUM/

-   `scripts/` --- Differential expression, clustering, and visualization analyses\
-   `gprofiler2_ORSUM/` --- Pathway enrichment and ORSUM summarization workflows

------------------------------------------------------------------------

# Reproducing Differential Expression and Clustering Analyses

This section describes how to reproduce:

-   Differential expression analyses (LRT and Wald tests by stage; DESeq2-based)
-   LRT clustering analysis using DEGreport

## 1. Clone the Repository

git clone https://github.com/cbib/cSCC_continuum_analyses\
cd cSCC_continuum_analyses

## 2. Create the Conda Environment

conda create --file geomx.yml\
conda activate geomx_env

## 3. Configure Input and Output Paths

Edit the following scripts:

-   `GeoMx_lrt_reproduce.R`
-   `HTGseq_lrt_reproduce.R`

Modify the variables:

* path_to_input\
* path_to_output

### Required Input Files

1.  Download the processed datasets from GEO (see *Data Availability*).
2.  Place all required input files in the directory specified by `path_to_input`.
3.  Download the file: Heatmaps_gene_list.xlsx and place it in the same input directory. ############################################################################

⚠️ The scripts will not run without these files.

## 4. Run the Analyses

Rscript GeoMx_lrt_reproduce.R\
Rscript HTGseq_lrt_reproduce.R

------------------------------------------------------------------------

# Output Structure

All results will be written to the directory specified in
`path_to_output` with the following structure:

* DE_genes_by_stage/\
* Heatmap_plots/\
* Single_gene_plots/

### Contents

-   **DE_genes_by_stage/**\
    DESeq2 pairwise comparisons for:
    -   HTG-seq\
    -   GeoMx Macrophages\
    -   GeoMx PanCK
-   **Heatmap_plots/**\
    Heatmaps summarizing pathway enrichment results obtained using manual curation of gprofiler2 and ORSUM output.
-   **Single_gene_plots/**\
    Strip plots for candidate stage-specific cSCC biomarkers identified in the study.

Additionally:

-   Likelihood Ratio Test (LRT) plots\
-   `.csv` files containing gene groups

are saved directly in `path_to_output`.

------------------------------------------------------------------------

# Reproducing gprofiler2 and ORSUM Pathway Analysis

⚠️ These scripts are configured to run on a **SLURM-based HPC cluster**.

If you wish to reproduce:

1.  gprofiler2 pathway enrichment\
2.  ORSUM summarization

follow the steps below.

## 1. Create the Environment

conda create --file gprofiler2_ORSUM/gprofiler2_orsum.yml\
conda activate gprof_orsum_env

## 2. Input Requirements

The script:

gprofiler2_ORSUM/gprof_orsum.sh

expects a `.tsv` file as input and calls:

gprofiler2_ORSUM/gprofiler_standalone.R

Notes:

-   LRT and DEGreport outputs are saved as `.csv`.
-   For GeoMx PanCK results, only a subset of gene groups was used for pathway analysis.
-   ORSUM requires `.gmt` files.

## 3. Download GMT Files

Download `.gmt` files from:

https://biit.cs.ut.ee/gprofiler/gost

Before running ORSUM:

-   Rename file prefixes by replacing `:` with `_`
    -   Example: `GO:BP` → `GO_BP`

## 4. Run on SLURM

cd gprofiler2_ORSUM/

sbatch gprof_orsum.sh -i Publication_HTG_Edgeseq_gene_groups.tsv -o
/results -gmt /gprofiler/gmt/hsapiens -os '500' -org 'hsapiens'

sbatch gprof_orsum.sh -i Publication_PanCK_gene_groups.tsv -o /results
-gmt /gprofiler/gmt/hsapiens -os '500' -org 'hsapiens'

sbatch gprof_orsum.sh -i Publication_Macrophages_gene_groups.tsv -o
/results -gmt /gprofiler/gmt/hsapiens -os '500' -org 'hsapiens'

Parameters:

-   `-i` : input gene group file\
-   `-o` : output directory\
-   `-gmt` : path to GMT files\
-   `-os` : ORSUM size threshold\
-   `-org` : organism (e.g., `hsapiens`)

------------------------------------------------------------------------

# Contact and Support

The code in this repository was developed by\
**Sergio Oterino-Sogo**

LinkedIn:\
https://www.linkedin.com/in/sergio-oterino-sogo-phd-181962164/

For reproducibility issues, please open a GitHub issue:\
https://github.com/cbib/cSCC_continuum_analyses/issues
