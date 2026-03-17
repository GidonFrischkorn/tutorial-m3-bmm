# Fitting the Memory Measurement Model (M3) with the bmm R Package: A Tutorial

Companion repository for the tutorial paper on fitting M3 (Oberauer & Lewandowsky, 2019) using the [bmm](https://venpopov.github.io/bmm/) R package. The paper demonstrates four progressively complex applications — simple span, complex span, a custom model with separate filtering parameters, and parameter recovery simulation.

## Repository Structure

This GitHub repository contains the source code and data needed to reproduce all analyses. Rendered outputs (fitted models, figures, manuscript PDFs) are available on [OSF](https://osf.io/yb7wm/).

### Included in the GitHub repository

```text
├── manuscript/          Quarto manuscript source (.qmd) and references
│   ├── tutorial-m3-bmm.qmd
│   └── references.bib
│
├── scripts/             Standalone R scripts for each tutorial
│   ├── tutorial1_simple_span.R              Simple span M3 (ss)
│   ├── tutorial2_complex_span.R             Complex span M3 (cs)
│   ├── tutorial3_custom_filtering.R         Custom M3 — separate ra/rc
│   ├── tutorial4_parameter_recovery.R       Parameter recovery simulation
│   ├── tutorial3_parameter_recovery_simple.R Appendix — memory updating (simple)
│   ├── tutorial3_parameter_recovery.R       Appendix — memory updating (full)
│   ├── prepare_Li2026_data.R                Data preparation for Tutorial 2
│   ├── figure_task_example.R                Figure 1: task diagram
│   ├── figure_m3_activations.R              Figure 2: activation decomposition
│   └── figure_m3_schematic.R                M3 framework schematic
│
├── data/                Experimental datasets
│   ├── Oberauer_2019_SimpleSpan_Exp1.dat    Tutorial 1: trial-level data
│   ├── Oberauer_2019_SimpleSpan_Exp2.dat    Tutorial 1: trial-level data
│   ├── Oberauer_2019_SimpleSpan_agg.csv     Tutorial 1: aggregated
│   ├── Li_2026_ComplexSpan_Exp1.csv         Tutorial 2: trial-level data
│   └── Li_2026_ComplexSpan_Exp1_agg.csv     Tutorial 2: aggregated
│
└── functions/           Shared plotting utilities
    └── clean_plot.R
```

### Available on OSF only

The following directories are excluded from the GitHub repository because they contain large binary files. They are archived on [OSF](https://osf.io/yb7wm/):

```text
├── output/              Cached model fits (.rds) — fitted brms/bmm objects
├── figures/             Generated figures (.pdf) — all manuscript figures
└── manuscript/*.pdf     Rendered manuscript (PDF, DOCX)
```

To reproduce locally, run the scripts in `scripts/` — fitted models will be cached in `output/` via the `file` argument in `bmm()`.

## Tutorials

**Tutorial 1 — Simple Span** (`tutorial1_simple_span.R`): Introduces the complete M3 workflow using data from Oberauer (2019). Covers model specification with `m3(version = "ss")`, fitting with `bmm()`, posterior predictive checks, parameter interpretation, choice rule comparison (simple vs. softmax), and hypothesis testing with `brms::hypothesis()`.

**Tutorial 2 — Complex Span** (`tutorial2_complex_span.R`): Extends the workflow to a task with distractors using data from Li, Frischkorn, & Oberauer (2026). Demonstrates `m3(version = "cs")`, the distractor filtering parameter `f`, handling non-identified parameters via constant priors, and condition-level hypothesis tests.

**Tutorial 3 — Custom Filtering** (`tutorial3_custom_filtering.R`): Defines a custom M3 with separate ratio parameters for item memory (`ra`) and context binding (`rc`) filtering, using the same data as Tutorial 2. Demonstrates user-defined activation formulas, link functions for bounded parameters, and model comparison with the standard `cs` version via bridge sampling.

**Tutorial 4 — Parameter Recovery** (`tutorial4_parameter_recovery.R`): Simulates data with known parameters using `rm3()`, fits the custom filtering model, and evaluates parameter recovery across a grid of sample sizes and trial counts. Demonstrates how to assess model identifiability and plan experimental designs.

**Appendix — Memory Updating** (`tutorial3_parameter_recovery_simple.R` and `tutorial3_parameter_recovery.R`): Supplementary scripts defining a custom M3 for a memory updating task. The simplified script walks through a single simulation cell; the full script varies sample size and trials per condition across a 3 × 3 design grid.

## Getting Started

Install required packages:

```r
install.packages("pacman")
pacman::p_load(here, bmm, brms, cmdstanr, tidyverse, tidybayes, patchwork, gghalves)
```

[cmdstanr](https://mc-stan.org/cmdstanr/) requires a working CmdStan installation. See `cmdstanr::install_cmdstan()`.

Open the RStudio project (`tutorial-m3-bmm.Rproj`) and run any script from `scripts/`. Fitted models are cached in `output/` via the `file` argument in `bmm()`, so re-running a script will load cached fits rather than re-fitting.

## Data Sources

- **Oberauer (2019)**: Oberauer, K. (2019). Working memory capacity limits memory for bindings. *Journal of Cognition, 2*(1), 40. https://doi.org/10.5334/joc.86. Data: [https://osf.io/vekpd/](https://osf.io/vekpd/)
- **Li et al. (2026)**: Li, C., Frischkorn, G. T., & Oberauer, K. (2026). Can we process information without encoding it into working memory? *Journal of Experimental Psychology: Learning, Memory, and Cognition*. https://doi.org/10.1037/xlm0001585. Data: [https://osf.io/wpcx5/](https://osf.io/wpcx5/)

## OSF Repository

All data, code, and materials are also archived on OSF: [https://osf.io/yb7wm/](https://osf.io/yb7wm/)

## References

- Oberauer, K., & Lewandowsky, S. (2019). Simple measurement models for complex working-memory tasks. *Psychological Review, 126*(6), 880–932. https://doi.org/10.1037/rev0000159
- Oberauer, K. (2019). Working memory capacity limits memory for bindings. *Journal of Cognition, 2*(1), 40. https://doi.org/10.5334/joc.86
- Li, C., Frischkorn, G. T., & Oberauer, K. (2026). Can we process information without encoding it into working memory? *Journal of Experimental Psychology: Learning, Memory, and Cognition*. https://doi.org/10.1037/xlm0001585
- Frischkorn, G. T. & Popov, V. (2025). A tutorial for estimating Bayesian hierarchical mixture models for visual working memory tasks. *Behavior Research Methods*.
- bmm package: [https://venpopov.github.io/bmm/](https://venpopov.github.io/bmm/)
