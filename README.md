# 🧬 Heterozygosity Visualization in Brassica rapa

This repository contains R scripts and data used to visualize heterozygous-like genomic regions in a Brassica rapa sequencing dataset.  
The analysis identifies genomic segments that deviate from expected homozygosity using allele-balance (AB) information derived from read depths.

---

## 📂 Overview

When a sample is expected to be homozygous, most variant sites should show only one allele — either reference or alternative.  
To detect potential heterozygous or mixed regions (caused by structural variation, mapping bias, or residual heterozygosity), this pipeline:

1. Calculates allele balance (AB) for each variant:
   AB = altDP / (refDP + altDP)
2. Flags sites with 0.3 ≤ AB ≤ 0.7 as heterozygous-like.
3. Aggregates results in 250 kb windows across each chromosome.
4. Computes % het-like per window:
   % het-like = (number of het-like sites / total sites) × 100
5. Generates a chromosome-scale ideogram showing:
   - Centromere positions (green)
   - Chromosome arms (white)
   - Grey “alarm” regions (windows with ≥ 5% het-like sites)
   - Blue gradient track (right edge): % het-like per window

---

## 🧪 Input

A tab-delimited text file containing variant-level information. Minimum required columns:

| CHROM | POS | REF | ALT | refDP | altDP | QUAL | total (optional) |
|-------|-----|-----|-----|-------|-------|------|------------------|
| A01   | 3045| G   | C   | 6     | 6     | 85   | 12               |
| A01   | 3093| T   | C   | 0     | 20    | 225  | 20               |
| ...   | ... | ... | ... | ...   | ...   | ...  | ...              |

- CHROM: chromosome name (e.g., A01)
- POS: position in bp
- refDP: reads supporting reference allele
- altDP: reads supporting alternative allele
- QUAL: variant quality score
- total: total read depth (optional; if missing the script can compute it as refDP + altDP)

---

## ⚙️ Requirements

R (≥ 4.0) and the following R packages:

```r
install.packages(c(
  "data.table", "dplyr", "tidyr", "forcats",
  "ggplot2", "ggideogram", "ggnewscale", "grid", "scales"
))
```

(If any package is not on CRAN, install from the appropriate source, e.g., Bioconductor or GitHub.)

---

## ▶️ Usage

1. Put your input file (example: `filter_uniq_paired_turnip15.txt`) in the working directory.
2. Run the plotting script:

```bash
Rscript plot_heterozygosity_ggideogram.R
# or, if the script accepts an argument:
# Rscript plot_heterozygosity_ggideogram.R filter_uniq_paired_turnip15.txt
```

Outputs:
- `heterozygosity_turnip15_ggideo.pdf` — high-quality figure for publication
- `heterozygosity_turnip15_ggideo.png` — lightweight preview for sharing

---

## 📊 Output interpretation

Each vertical bar in the ideogram represents one A-genome chromosome (A01–A10).

Legend:
- 🟩 Green — Centromere region (acen)
- ⚪ White — Chromosome arms (gneg)
- 🔵 Blue gradient (right strip) — % of heterozygous-like sites per 250 kb window
- ⚫ Grey shading — Alarm window (≥ 5% het-like sites)
- • Small dots (left edge) — Individual heterozygous-like sites

Axes:
- Y-axis: Chromosome position in Mb (centromere aligned to 0)
- X-axis: Chromosomes A01–A10

---

## 🧠 Biological interpretation

- Low % het-like (white / very light blue): Region is homozygous as expected.
- High % het-like (dark blue / grey): Region shows potential heterozygosity — possible causes:
  - True residual heterozygosity
  - Misalignment or duplicated regions
  - Structural variation or CNV
  - Contamination or mixed samples

Follow-up analyses:
- Inspect raw read alignments (IGV)
- Generate per-region read depth plots
- Run dedicated CNV or structural variant callers

---

## ⚙️ Parameters you may want to change

These are typically defined at the top of the R script:
- het_low / het_high: heterozygous-like AB interval (default 0.3–0.7)
- alarm_pct_het: alarm threshold for a window (default 5)
- min_QUAL: minimum variant quality to keep
- min_DP: minimum total depth to keep
- window_size: window size for aggregation (default 250000)

---

## 🪄 Notes & tips

- The script automatically handles windowing and centromere alignment for the A-genome chromosomes.
- If `total` is missing, it is computed as refDP + altDP.
- Adjust the AB thresholds and alarm percentile to match your experiment and expected error models.
- Filter low-quality or low-depth variants to reduce noise.

---

## 🧰 License

This project is released under the MIT License. 

