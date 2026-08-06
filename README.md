# ALPACA Pipeline

Nextflow DSL2 pipeline for tumor copy-number evolution analysis. Integrates
Purple CNV calls, [RefPhase](https://github.com/amcrabtree/refphase)
allele-specific phasing, and the [ALPACA](https://pypi.org/project/alpaca/)
optimizer to infer clonal copy-number profiles along a phylogenetic tree.

## Repository layout

```
ALPACA_pipeline.nf          ← main workflow
nextflow.config              ← default params, executor, conda, Gurobi env
local.config.example         ← template for your own environment (copy to local.config)
modules/                     ← one process per file
├── name-conversion.nf       ← NAME_CONVERSION
├── harmonize.nf             ← HARMONIZE
├── clonal-info.nf           ← CLONAL_INFO
├── refphase.nf               ← REFPHASE
├── calculate-ci.nf           ← CALCULATE_CI
├── alpaca.nf                 ← ALPACA (core solver)
├── get-stats.nf              ← GET_STATS (R diagnostic plots)
└── post-process.nf           ← POST_PROCESS (CCD + WGD + plot-tumour)
assets/                       ← R helper scripts sourced by modules
old_pipelines/                ← superseded monolithic versions
```

## Pipeline steps

```
sample_info TSV + batch_info TSV + Purple output + SNP BAF files
        │
        ▼
NAME_CONVERSION   → {pat}_cf_converted.tsv, _segments_converted.tsv, _purity_ploidy_converted.tsv
        │
        ├──────────────────────┐
        ▼                      ▼
HARMONIZE                 CLONAL_INFO
→ {pat}-segments.tsv      → tree_paths.json
→ {pat}-snps.tsv          → cp_table.csv
→ {pat}-purity_ploidy.tsv
        │
        ▼
REFPHASE
→ {pat}-refphase-segmentation.tsv
→ {pat}-refphase-phased-snps.tsv.gz
→ {pat}-refphase-sample-data-updated.tsv
        │
        ▼
CALCULATE_CI
→ ci_table.csv
→ ALPACA_input_table.csv
        │
        ▼ (joined with CLONAL_INFO output)
ALPACA  [Gurobi, SLURM]
→ ALPACA_output_{pat}.csv
→ cn_change_to_ancestor.csv
→ *_report.csv, run_gap_summary.csv  (diagnostic, optional)
        │
        ├──────────────────────────────────────┐
        ▼                                      ▼
GET_STATS                                 POST_PROCESS
→ histograms/*.png                        → clone_copy_number_diversity_scores.csv
                                          → wgd_ratio_scores.csv
                                          → {pat}_plots.ipynb
```

## Requirements

| Tool | Notes |
|---|---|
| [Nextflow](https://www.nextflow.io/) | DSL2 |
| [ALPACA](https://pypi.org/project/alpaca/) Python package | provides the `alpaca` CLI (`run`, `ccd`, `wgd`, `plot-tumour`) |
| [Gurobi](https://www.gurobi.com/) optimizer | needs a valid license |
| [RefPhase](https://github.com/amcrabtree/refphase) | R package |
| conda/mamba | used to provision the environment via Nextflow's `conda` directive |
| SLURM | or adapt `nextflow.config` for another executor |

## Setup

1. Create a conda environment with `alpaca`, `refphase`, and their R/Python
   dependencies installed.
2. Install Gurobi and obtain a license.
3. Copy `local.config.example` to `local.config` and fill in the paths for
   your own environment (input data locations, conda env, Gurobi install,
   SLURM queue). `nextflow.config` automatically includes `local.config` if
   it exists, so there's no extra flag needed at run time. `local.config` is
   gitignored so your local paths never get committed.

## Inputs

| Parameter | Description |
|---|---|
| `--seg` | Purple segmentation TSV |
| `--pp` | Purple purity/ploidy estimates TSV |
| `--snp_path` | directory of per-sample SNP BAF files (`{patient}/{sample}.amber.baf.tsv[.gz]`) |
| `--batch_info` | TSV with columns `patient`, `batch`, `model`, `path_to_trees`, `file_cf`, `file_tree` |
| `--conversion_table` | sample-name conversion table (CSV) |
| `--sample_info` | TSV with `sample` and `patient` columns |
| `--pubDir` | output directory |

## Run command

```bash
nextflow run ALPACA_pipeline.nf \
    -resume \
    -profile conda \
    --sample_info /path/to/sample_info.tsv \
    --pubDir /path/to/output
```

## Output per patient (`{pubDir}/{patient}/`)

| File | Produced by |
|---|---|
| `ALPACA_output_{pat}.csv` | ALPACA — predicted CN per clone per segment |
| `cn_change_to_ancestor.csv` | ALPACA — CN changes along tree edges |
| `ci_table.csv` | CALCULATE_CI — confidence intervals per segment |
| `ALPACA_input_table.csv` | CALCULATE_CI — reformatted CN table for ALPACA |
| `cp_table.csv` | CLONAL_INFO — clone proportions per sample |
| `tree_paths.json` | CLONAL_INFO — phylogenetic tree paths |
| `{pat}-segments.tsv` | HARMONIZE — unified CN segments |
| `{pat}-snps.tsv` | HARMONIZE — BAF/logR SNP data |
| `{pat}-purity_ploidy.tsv` | HARMONIZE — purity/ploidy per sample |
| `{pat}-refphase-*.tsv(.gz)` | REFPHASE — phased allele-specific CN |
| `clone_copy_number_diversity_scores.csv` | POST_PROCESS (`alpaca ccd`) |
| `wgd_ratio_scores.csv` | POST_PROCESS (`alpaca wgd`) |
| `{pat}_plots.ipynb` | POST_PROCESS (`alpaca plot-tumour --plot_output_mode notebook`) |
| `*_report.csv`, `run_gap_summary.csv` | ALPACA — diagnostic reports |
| `histograms/*.png` | GET_STATS — RefPhase vs ALPACA comparison plots |

## Known quirks

- **Chromosome X in segment names**: `alpaca ccd` and `alpaca wgd` require
  purely numeric chromosome identifiers (`^\d+_\d+_\d+$`). The ALPACA output
  contains `X_start_end` segments, which fail validation. `POST_PROCESS`
  rewrites them to `23_start_end` in a temporary file
  (`alpaca_output_recoded.csv`) before passing to those two commands.
  `alpaca plot-tumour` receives the original file with `X_` names.

- **PDF plots require Chrome**: `alpaca plot-tumour --plot_output_mode pdf`
  uses Kaleido, which requires Google Chrome. If Chrome isn't available in
  your environment, generate plots as Jupyter notebooks
  (`--plot_output_mode notebook`) instead and open the `.ipynb` output
  locally to view interactive figures.
