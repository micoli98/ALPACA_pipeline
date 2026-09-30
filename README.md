# ALPACA-SVclone Pipeline

Nextflow DSL2 pipeline for tumor copy-number evolution analysis. Integrates
Purple CNV calls, [RefPhase](https://github.com/amcrabtree/refphase)
allele-specific phasing, and the [ALPACA](https://pypi.org/project/alpaca/)
optimizer to infer clonal copy-number profiles along a phylogenetic tree, and
runs [SVclone](https://github.com/mcmero/SVclone) to cluster structural
variants by cancer cell fraction.

## Repository layout

```
main.nf                        ← top-level entrypoint: runs ALPACA_SUBWORKFLOW + SVCLONE_SUBWORKFLOW, joins output by patient
nextflow.config                ← default params, executor, conda, Gurobi env
local.config.example           ← template for your own environment (copy to local.config)
subworkflows/
├── alpaca.nf                  ← ALPACA_SUBWORKFLOW: clonal copy-number evolution (per patient)
└── svclone.nf                 ← SVCLONE_SUBWORKFLOW: SV clustering by cancer cell fraction (per sample)
modules/                       ← one process per directory (modules/<ProcessName>/main.nf)
├── NameConversion/            ← NAME_CONVERSION
├── Harmonize/                 ← HARMONIZE
├── ClonalInfo/                ← CLONAL_INFO
├── Refphase/                  ← REFPHASE
├── CalculateCi/               ← CALCULATE_CI
├── Alpaca/                    ← ALPACA (core solver)
├── GetStats/                  ← GET_STATS (R diagnostic plots)
├── PostProcess/               ← POST_PROCESS (CCD + WGD + plot-tumour, v0.3.1)
├── GetPatientInfo/            ← GET_PATIENT_INFO (currently unused/unwired)
├── GetInsertStd/              ← GET_INSERT_STD (SVclone)
├── Config/                    ← CONFIG (SVclone)
├── FilterSV/                  ← FILTER_SV (SVclone)
├── Annotation/                ← ANNOTATION (SVclone)
├── Count/                     ← COUNT (SVclone)
├── PrepareFilterInput/        ← PREPARE_FILTER_INPUT (SVclone)
├── Filter/                    ← FILTER (SVclone)
└── Clustering/                ← CLUSTERING (SVclone)
assets/                        ← R helper scripts sourced by ALPACA modules + convert_refphase.py (CALCULATE_CI)
old_pipelines/                 ← superseded monolithic versions
```

## Pipeline steps

`main.nf` runs both subworkflows over the same cohort (`sample_info`) and
joins their outputs by patient into one terminal in-memory channel. There is
no combined `publishDir` for that join — each process publishes its own
files under `{outdir}/{patient}/alpaca/` (ALPACA side) or
`{outdir}/{patient}/svclone/{sample}/` (SVclone side) as it runs.

### ALPACA subworkflow (per patient)

```
sample_info TSV + batch_info TSV + Purple output + AMBER BAF + COBALT ratio files
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
REFPHASE  (cn_method = "none", phase_all_segs = FALSE)
→ {pat}-refphase-segmentation.tsv
→ {pat}-refphase-phased-snps.tsv.gz
→ {pat}-refphase-sample-data-updated.tsv
        │
        ▼
CALCULATE_CI  (bootstrap CIs, assets/convert_refphase.py; needs cp_table.csv)
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

Notes on the copy-number inputs:

- **BAF** — `HARMONIZE` takes the raw AMBER `tumorBAF` column (not
  `tumorModifiedBAF`, which AMBER folds to ≥ 0.5 and so carries no phase
  information). Together with `cn_method = "none"` this fixed the ploidy
  inflation seen in earlier runs.
- **logR** — GC-corrected logR from COBALT:
  `log2(tumorGCRatio / referenceGCDiploidRatio)` of the 1 kb COBALT window
  containing each SNP. SNPs in masked (`-1`) or zero-coverage windows get no
  logR and are dropped (RefPhase rejects NA logR); the percentage dropped per
  sample is printed in the task log.
- **RefPhase** runs with `cn_method = "none"`, so Purple's allele-specific
  copy numbers are phased but not re-estimated.
- **Confidence intervals** — `CALCULATE_CI` runs ALPACA's RefPhase converter
  (`assets/convert_refphase.py`, a wrapper around ALPACA's
  `convert_refphase_output` functions): for each segment × sample, SNPs are
  resampled 1000 times (`--n_bootstrap`) and the 95% CI of the ASCAT copy
  number from phased BAF + logR is taken, centred on the RefPhase copy
  number. `--heterozygous_SNPs_threshold 0` keeps every segment; segments
  with ≤ 10 het SNPs, no SNPs, or no usable logR get a fixed CI span of 0.5.
  Only samples present in CONIPHER's `cp_table.csv` are kept.

### SVclone subworkflow (per sample)

```
sample_info TSV + batch_info TSV + coverage_info TSV + Purple output (seg/pp/SV VCF) + BAM/BAI
        │
        ├──────────────────────────────┐
        ▼                              ▼
GET_INSERT_STD                    (purity filter: only samples
→ insert size stdev                above --min_purity proceed)
        │                              │
        └──────────────┬───────────────┘
                        ▼
                     CONFIG
                     → svclone_config.ini
        │
        ▼
FILTER_SV
→ {sample}.purple.breakpoints.vcf
        │
        ▼
ANNOTATION
→ {sample}_svin.txt, read_params.txt
        │
        ▼
COUNT
→ {sample}_svinfo.txt
        │
        ▼ (joined with PREPARE_FILTER_INPUT output)
PREPARE_FILTER_INPUT
→ {sample}_snvs_for_svclone.vcf, {sample}_ascat.csv, {sample}_pp.tsv
        │
        ▼
FILTER
→ {sample}_filtered_svs.tsv, {sample}_filtered_snvs.tsv, purity_ploidy.txt
        │
        ▼
CLUSTERING
→ SVclone cluster output files (per sample)
```

## Requirements

| Tool | Notes |
|---|---|
| [Nextflow](https://www.nextflow.io/) | DSL2 |
| [ALPACA](https://pypi.org/project/alpaca/) Python package | provides the `alpaca` CLI (`run`, `ccd`, `wgd`, `plot-tumour`); academic/non-commercial license |
| [Gurobi](https://www.gurobi.com/) optimizer | needs a valid license |
| [RefPhase](https://bitbucket.org/schwarzlab/refphase) | R package, not on CRAN/conda — installed via `devtools::install_bitbucket()` |
| [SVclone](https://github.com/mcmero/SVclone) | provides the `svclone` CLI (`annotate`, `count`, `filter`, `cluster`); older Python 3.6 environment |
| conda/mamba | used to provision the environments via Nextflow's `conda` directive |
| SLURM | or adapt `nextflow.config` for another executor |

## Setup

1. Create the conda environment from `alpaca_environment.yml`:

   ```bash
   conda env create -f alpaca_environment.yml
   ```

   The `ANNOTATION`, `COUNT`, `PREPARE_FILTER_INPUT`, `FILTER` and `CLUSTERING` processes (SVclone
   subworkflow) run in a separate, older (Python 3.6) environment instead — create it from
   `svclone_environment.yml`:

   ```bash
   conda env create -f svclone_environment.yml
   ```

2. Install RefPhase into that environment (not distributed via conda/CRAN):

   ```bash
   conda run -n alpaca Rscript -e 'devtools::install_bitbucket("schwarzlab/refphase")'
   ```

3. Install Gurobi and obtain a license (the `gurobi` conda package installs
   the Python bindings; you still need a license file, see
   [gurobi.com](https://www.gurobi.com/downloads/)).
4. Copy `local.config.example` to `local.config` and fill in the paths for
   your own environment (input data locations, conda envs, Gurobi install,
   SLURM queue). `nextflow.config` automatically includes `local.config` if
   it exists, so there's no extra flag needed at run time. `local.config` is
   gitignored so your local paths never get committed.

## Inputs

| Parameter | Description |
|---|---|
| `--seg` | Purple segmentation TSV |
| `--pp` | Purple purity/ploidy estimates TSV |
| `--snp_path` | directory of per-sample AMBER BAF and COBALT ratio files (`{patient}/{sample}.amber.baf.tsv[.gz]`, `{patient}/{sample}.cobalt.ratio.tsv.gz`) |
| `--batch_info` | TSV with columns `patient`, `batch`, `model`, `path_to_trees`, `file_cf`, `file_tree`, `snv_vcf` |
| `--conversion_table` | sample-name conversion table (CSV) |
| `--sample_info` | TSV with `sample` and `patient` columns |
| `--outdir` | output directory |
| `--bam_dir` | directory of per-sample BAM/BAI files (SVclone) |
| `--purple_dir` | directory of per-patient Purple SV VCFs (SVclone) |
| `--coverage_info` | TSV with columns `sample`, `meanCoverage`, `MEAN_READ_LENGTH`, `MEAN_INSERT_SIZE` (SVclone) |
| `--conversion_script` | R script sourced by `PREPARE_FILTER_INPUT` for mutation-key matching (SVclone) |
| `--min_purity` | Purple purity threshold below which samples are excluded from SVclone (default `0.1`) |

## Run command

```bash
nextflow run main.nf \
    -resume \
    -profile conda \
    --sample_info /path/to/sample_info.tsv \
    --outdir /path/to/output
```

## Output

### Per patient (`{outdir}/{patient}/alpaca/`) — ALPACA subworkflow

| File | Produced by |
|---|---|
| `ALPACA_output_{pat}.csv` | ALPACA — predicted CN per clone per segment |
| `cn_change_to_ancestor.csv` | ALPACA — CN changes along tree edges |
| `ci_table.csv` | CALCULATE_CI — bootstrap confidence intervals per segment and sample |
| `ALPACA_input_table.csv` | CALCULATE_CI — reformatted CN table for ALPACA |
| `cp_table.csv` | CLONAL_INFO — clone proportions per sample |
| `tree_paths.json` | CLONAL_INFO — phylogenetic tree paths |
| `{pat}-segments.tsv` | HARMONIZE — unified CN segments |
| `{pat}-snps.tsv` | HARMONIZE — AMBER BAF + COBALT logR per SNP |
| `{pat}-purity_ploidy.tsv` | HARMONIZE — purity/ploidy per sample |
| `{pat}-refphase-*.tsv(.gz)` | REFPHASE — phased allele-specific CN |
| `clone_copy_number_diversity_scores.csv` | POST_PROCESS (`alpaca ccd`) |
| `wgd_ratio_scores.csv` | POST_PROCESS (`alpaca wgd`) |
| `{pat}_plots.ipynb` | POST_PROCESS (`alpaca plot-tumour --plot_output_mode notebook`) |
| `*_report.csv`, `run_gap_summary.csv` | ALPACA — diagnostic reports |
| `histograms/*.png` | GET_STATS — RefPhase vs ALPACA comparison plots |

### Per sample (`{outdir}/{patient}/svclone/{sample}/`) — SVclone subworkflow

| File | Produced by |
|---|---|
| `svclone_config.ini` | CONFIG — per-sample SVclone configuration |
| `{sample}.purple.breakpoints.vcf` | FILTER_SV — breakpoint-only SV VCF |
| `{sample}_svin.txt`, `read_params.txt` | ANNOTATION — annotated SV input / read parameters |
| `{sample}_svinfo.txt` | COUNT — read-support counts per SV |
| `{sample}_snvs_for_svclone.vcf`, `{sample}_ascat.csv`, `{sample}_pp.tsv` | PREPARE_FILTER_INPUT — reformatted SNV/CNV/purity-ploidy input |
| `{sample}_filtered_svs.tsv`, `{sample}_filtered_snvs.tsv`, `purity_ploidy.txt` | FILTER — filtered SVs/SNVs |
| SVclone cluster output | CLUSTERING (`svclone cluster`) — CCF clusters |

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

- **Low-purity samples are skipped in SVclone**: `SVCLONE_SUBWORKFLOW` reads
  Purple purity from `--pp` and drops any sample at or below `--min_purity`
  (default `0.1`) before `GET_INSERT_STD` runs, so those samples produce no
  SVclone output.

- **`CALCULATE_CI` is slow**: ALPACA's bootstrap is a pandas loop (1000
  resamples per segment × sample, single core). Whole-genome patients with
  ~16k segment×sample groups take 1–15 h per patient; the process is given
  `time = '48h'` and `memory = 32.GB` in `nextflow.config`.
