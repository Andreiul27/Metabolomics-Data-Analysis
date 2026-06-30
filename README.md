# Untargeted Metabolomics Group Comparison Analysis

An R Markdown pipeline for analyzing an untargeted LC-MS metabolomics dataset across three experimental/treatment groups, with paired statistical testing, volcano plots, exploratory clustering, and PCA (with and without pooled QC samples).

The script currently runs on a **simulated dataset** (generated in-script) so it can be shared and tested without exposing real patient data. 
## Workflow

The analysis, in order, covers:

1. **Import** — reads an `.xlsx` feature table exported from LC-MS processing software.
2. **Reshape** — converts the wide feature table (one row per feature, one column per sample) into long format, parsing sample column names into `Patient` and `Group` identifiers.
3. **Duplicate check** — a small utility function to confirm no duplicate rows exist after reshaping.
4. **Statistical testing** — paired t-tests (Group 1 vs 2, 1 vs 3, 2 vs 3) with Benjamini–Hochberg FDR correction, computed three ways:
   - on group-level means per feature
   - per patient, across all features
   - per feature, including log2 fold change and -log10(p-value) for downstream volcano plotting
5. **Volcano plot** — interactive (Plotly) scatter plot of log2 fold change vs. -log10(p-value), flagging features as significant at p ≤ 0.05 and |log2FC| > 1.
6. **Boxplots / violin plots** — per-feature distribution plots for significant features (disabled by default check [Notes](#notes)).
7. **Unsupervised clustering** — k-means clustering of features across patients, with and without a high-intensity outlier feature excluded.
8. **PCA** — principal component analysis of samples, run both with and without the pooled QC (PQC) samples included, visualized with `factoextra`/Plotly.

## Requirements

R (≥ 4.2 recommended) with the following packages:

```r
install.packages(c(
  "readxl", "reshape2", "dplyr", "ggplot2", "ggfortify", "tidyr",
  "purrr", "broom", "RColorBrewer", "plotly", "factoextra",
  "corrplot", "tidyverse", "viridis"
))
```

## Data input

The script expects an Excel feature table with:

- `ID` — unique feature identifier
- `m/z`, `RT [min]` — mass-to-charge ratio and retention time
- `Blank` — blank sample intensity
- Sample columns named `<patient_number>_<group_number>` (e.g. `1_1`, `1_2`, `1_3`, `2_1`, ...) for 20 patients × 3 groups
- QC columns named `QC_<n>` for pooled QC injections



Then comment out or remove the "Simulate the data" chunk to use your real import instead of the synthetic dataset.

## Usage

1. Open `metabolomics_group_comparison_analysis.Rmd` in RStudio.
2. Update the `path` variable (and re-enable the real `read_excel` import) if running on real data.
3. Knit to HTML, or run chunks interactively in the console.

## Notes

- Several chunks are intentionally commented out:
  - The per-feature boxplot/violin loop (saves individual HTML widgets per significant feature)  you can uncomment to regenerate.
  - The per-feature m/z and RT scatter plot loops are **not** meant to be run as is; they can take several hours given the full feature set. Re-enable only if needed, and consider running on a subset first.
- `num_clusters` controls the number of k-means clusters and can be adjusted.
- The PCA-with-QC section recenters PC1/PC2 using fixed offset values derived from a prior run on real data; revisit these if your dataset's variance structure differs substantially.

