# dNdScv_R_analysis

**VCF → QC tables → Publication-grade plots → Shiny dashboard →
dNdScv-ready mutation table (R + bcftools-backed)**

`dNdScv_R_analysis` is an R-based downstream analysis module that:

-   Reads one or multiple VCFs via **bcftools query** (robust and
    lightweight; avoids heavy VCF parsing dependencies)
-   Produces a **dNdScv-ready mutation table**: `dndscv_input.csv`
-   Generates **quality control summaries and publication-grade plots**
-   Optionally scaffolds a **Shiny dashboard** for interactive
    inspection

This repository demonstrates:

-   Reproducible R environments using `renv`
-   Efficient data processing pipelines
-   Structured CLI-driven workflows
-   Publication-ready visualisation with `ggplot2`
-   Interactive reporting with Shiny

------------------------------------------------------------------------

## Relationship to dNdScv_prep (Python upstream workflow)

This repository complements the upstream Python pipeline:

**dNdScv_prep**\
FASTQ → BAM → VCF → mutation tables\
GitHub: https://github.com/samson-olofinsae/dNdScv_prep

### Recommended toolchain

1.  Generate VCFs using `dNdScv_prep` (or bring externally generated
    VCFs)
2.  Run `dNdScv_R_analysis` to produce QC, visualisations, and
    `dndscv_input.csv`
3.  Run the Martincorena lab `dndscv` package downstream

------------------------------------------------------------------------

## Requirements

### System

-   **R ≥ 4.2** (tested on R 4.4.0)
-   **bcftools ≥ 1.10** (tested on 1.21)

Check versions:

``` bash
bcftools --version | head
R --version | head
```

### R packages (managed via renv)

-   optparse
-   data.table / dplyr
-   ggplot2
-   shiny

This repository includes `renv.lock` to ensure environment
reproducibility.

------------------------------------------------------------------------

## Quick Start

Clone:

``` bash
git clone https://github.com/samson-olofinsae/dNdScv_R_analysis.git
cd dNdScv_R_analysis
```

Restore environment (first run only):

``` bash
R -q -e "install.packages('renv'); renv::restore(prompt=FALSE)"
```

------------------------------------------------------------------------

## Preparing Your VCF Folder

Copy VCF files into:

``` bash
data/user_vcf/
```

Example:

``` bash
cp /path/to/your/*.vcf.gz data/user_vcf/
cp /path/to/your/*.vcf.gz.tbi data/user_vcf/  # optional but recommended
```

Notes:

-   Input files must be `.vcf.gz` (bgzipped)
-   Index files (`.tbi`) are recommended

------------------------------------------------------------------------

## Running the Pipeline

``` bash
Rscript dndscv_R_analysis.R --vcf-glob "data/user_vcf/*.vcf.gz" --outdir dndscv_R_prep_output --write-qc --write-plots --write-shiny
```

Output structure:

    dndscv_R_prep_output/
    ├── dndscv_input.csv
    ├── qc_summary.csv
    ├── plots/
    │   ├── qc_variants_bar.png
    │   ├── qc_snv_indel_stacked.png
    │   ├── qc_snv_fraction.png
    │   └── qc_variants_by_chr.png
    └── shiny_app/
        └── app.R

Launch the Shiny dashboard:

``` bash
R -q -e "shiny::runApp('dndscv_R_prep_output/shiny_app')"
```

------------------------------------------------------------------------

## Publication-Grade Plots

The `plots/` directory contains publication-ready PNG figures suitable
for manuscripts, presentations, and panel review.

### qc_variants_bar.png

Total variant count per sample.\
Useful for identifying outlier samples or sequencing depth
inconsistencies.

### qc_snv_indel_stacked.png

Stacked SNV vs INDEL counts per sample.\
Highlights mutation spectrum balance and potential calling biases.

### qc_snv_fraction.png

Proportion of SNVs relative to total variants.\
Useful for quality assessment and mutation type profiling.

### qc_variants_by_chr.png

Distribution of variants across chromosomes.\
Helps detect chromosomal clustering, panel coverage bias, or artefacts.

------------------------------------------------------------------------

## plots_for_panel

To support panel review or interview demonstration:


plots_for_panel/

This folder can contain:

-   Real-world QC outputs from genuine datasets (anonymised)

- All data is anonymised and compliant with data governance policies.

------------------------------------------------------------------------

## Running dNdScv Downstream

Once `dndscv_input.csv` is generated:

``` r
library(dndscv)

m <- read.csv("dndscv_R_prep_output/dndscv_input.csv", stringsAsFactors = FALSE)
out <- dndscv(m)
out$sel_cv
```

------------------------------------------------------------------------

## Reproducibility

This repository uses **renv**:

-   `renv.lock` pins exact package versions
-   `renv::restore()` recreates the environment on any machine

------------------------------------------------------------------------

## License

MIT License
