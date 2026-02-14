# dNdScv_R_analysis
**VCF → QC tables/plots → Shiny dashboard → dNdScv-ready mutation table (R + bcftools-backed)**

`dNdScv_R_analysis` is an R-based downstream analysis module that:

- reads one or many VCFs via **bcftools query** (robust; avoids VariantAnnotation IO edge cases)
- produces a **dNdScv-ready** mutation table: `dndscv_input.csv`
- generates **QC summaries + plots**
- optionally writes a **Shiny** app scaffold for interactive inspection

This repo is designed to demonstrate R proficiency in:

- reproducible analysis structure (`renv`)
- `data.table` pipelines
- `ggplot2` visualisation
- Shiny interactive reporting

---

## Designed to complement: dNdScv_prep (Python)

This repo complements the upstream Python workflow:

**dNdScv_prep** (FASTQ → BAM → VCF → mutation tables)  
GitHub: https://github.com/samson-olofinsae/dNdScv_prep

**Recommended toolchain**

1) Use `dNdScv_prep` to generate VCFs (or bring your own VCFs)  
2) Use `dNdScv_R_analysis` to generate QC + visualisation + `dndscv_input.csv`  
3) Run dNdScv (Martincorena lab package) on the resulting table

---

## Requirements

### System
- **R** (>= 4.2 recommended; tested on **R 4.4.0**)
- **bcftools** in PATH (>= 1.10; tested on **1.21**)

Check:
```bash
bcftools --version | head
R --version | head
```

### R packages (managed by renv)
- optparse
- data.table
- ggplot2
- shiny

This repo includes `renv.lock` to reproduce the exact R package environment.

---

## Quick start

Clone:
```bash
git clone https://github.com/samson-olofinsae/dNdScv_R_analysis.git
cd dNdScv_R_analysis
```

Restore the R environment (first time only):
```bash
R -q -e "install.packages('renv'); renv::restore(prompt=FALSE)"
```

Run the pipeline on VCFs:
```bash
Rscript dndscv_R_analysis.R \
  --vcf-glob "data/vcfs/*.vcf.gz" \
  --outdir dndscv_R_prep_output \
  --write-qc \
  --write-plots \
  --write-shiny
```

Outputs:
```text
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
```

Run the Shiny app:
```bash
R -q -e "shiny::runApp('dndscv_R_prep_output/shiny_app')"
```

---

## Viewing the PNG plots

On Linux/WSL:
```bash
ls dndscv_R_prep_output/plots
xdg-open dndscv_R_prep_output/plots/qc_variants_bar.png
```

Open in Windows Explorer (from WSL):
```bash
explorer.exe dndscv_R_prep_output/plots
```

---

## dNdScv usage (downstream)

Once `dndscv_input.csv` is generated, you can run dNdScv:

```r
library(dndscv)
m <- read.csv("dndscv_R_prep_output/dndscv_input.csv", stringsAsFactors = FALSE)
out <- dndscv(m)
out$sel_cv
```

---

## Reproducibility

This repository uses **renv**:
- `renv.lock` records exact package versions
- `renv::restore()` recreates the environment on any machine

---

## License
MIT License.
