# mousetohuman_predict

Apply the mouse-to-human translatability models to your own murine blood
transcriptomic dataset.

This repository is a **standalone companion** to
[`animals_vax_atlas`](https://github.com/wapsyed/animals_vax_atlas). It contains
only what is needed to run the prediction notebook:

- `Run your analysis here/RunYourAnalysis_MouseToHuman.Rmd` — the step-by-step notebook.
- `Run your analysis here/example_dge_result.rds` — an example mouse DGE table (*S. aureus*).
- `scripts_notebooks/required.R` — minimal helpers (`theme_vaxgo`, `autoGSEA`).
- `tables/` — the annotation layers used as model features.
- `Modelling/Models/` — the trained random-forest models (best per task).

## Run it

### Option 1 — Binder (no installation)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/wapsyed/mousetohuman_predict/HEAD?urlpath=rstudio)

Click the badge, wait for the environment to build (the first launch takes a
few minutes), then open `Run your analysis here/RunYourAnalysis_MouseToHuman.Rmd`
and run the chunks. The environment is defined by `runtime.txt`, `install.R` and
`apt.txt`.

### Option 2 — Run locally

```r
# 1. Clone the repository
#    git clone https://github.com/wapsyed/mousetohuman_predict.git
# 2. Install the packages
source("install.R")
# 3. Open and knit the notebook
#    Run your analysis here/RunYourAnalysis_MouseToHuman.Rmd
```

## Input

Your input table must contain these columns:

`pathogen, treatment, timepoint, hgnc_symbol, mean_l2fc, ci_lower, ci_upper,
ave_expr, t, p_value, adj_p_val, b, se, sd, inverse_se, organism, condition`

The notebook runs GSEA against the BTM gene sets, builds the model features and
predicts the human response.

## Output — the four scores

Each gene is scored by a model trained without its condition (out-of-fold):

| Score | Meaning |
|:------|:--------|
| `score_shared` | Probability of being a shared leading-edge gene (mouse + human). |
| `score_rank` | Predicted human absolute rank (0–100). |
| `score_direction` | Probability of concordant direction (same sign in both species). |
| `score_translational` | `score_rank × score_direction`. |

## Models

The models were trained on mouse-to-human transfer across infection and injury
conditions using leave-one-pathogen-out cross-validation. A random forest was
the best algorithm for every task.

| File | Task |
|:-----|:-----|
| `Modelling/Models/rf_model_shared.rds` | Shared vs Mouse-only LEG classification |
| `Modelling/Models/rf_model_rank.rds` | Human absolute rank regression |
| `Modelling/Models/rf_model_direction.rds` | Directional concordance classification |

## Citation

> Prates-Syed WA, Lira AA, Cortes N, Silva JDQ, Hamaguchi B, Carvalho E,
> Castillo-Chávez A, Durães-Carvalho R, Cabral-Marques O, Sabino EC, Krieger JE,
> Hagan T, Cabral-Miranda G. *From Mice to Humans: Functional Modules Improve the
> Translatability of Transcriptomic Responses.* Genes and Immunity (under review).
