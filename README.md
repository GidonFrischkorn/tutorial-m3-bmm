# Applying the Memory Measurement Model (M3): A Tutorial Using the bmm R Package

Companion repository for the tutorial on fitting the Memory Measurement Model (M3; Oberauer & Lewandowsky, 2019) with the [bmm](https://venpopov.github.io/bmm/) R package. The tutorial works through four applications: a simple span task, a complex span task, a custom M3 with separate filtering parameters, and a parameter recovery simulation. Supplement 2 adds a parameter recovery study for a custom M3 of a memory updating task.

This repository holds the code, data, and manuscript sources. Fitted models and figures are archived on [OSF](https://osf.io/yb7wm/).

## Repository structure

```text
├── manuscript/          Quarto sources and references
│   ├── tutorial-m3-bmm.qmd                     Main text
│   ├── supplement1-methods.qmd                 Supplement 1: prior predictive checks, Bayes factors,
│   │                                           softmax as SDT, validating custom M3
│   ├── supplement2-parameter-recovery.qmd      Supplement 2: recovery for a memory updating M3
│   └── references.bib
│
├── scripts/
│   ├── prepare_Li2026_data.R                   Data preparation for Tutorials 2 and 3
│   ├── tutorial1_simple_span.R                 Tutorial 1: simple span (ss)
│   ├── tutorial2_complex_span.R                Tutorial 2: complex span (cs)
│   ├── tutorial3_custom_filtering.R            Tutorial 3: custom M3 with separate filtering (ra, rc)
│   ├── tutorial4_parameter_recovery.R          Tutorial 4: parameter recovery
│   ├── supplement1_loo_choice_rule.R           Supplement 1: LOO comparison of the choice rules
│   ├── supplement1_savage_dickey_stability.R   Supplement 1: Monte Carlo error of Savage-Dickey ratios
│   ├── supplement2_parameter_recovery_updating_simple.R   Supplement 2: single design cell
│   ├── supplement2_parameter_recovery_updating.R          Supplement 2: 3 x 3 design grid
│   ├── figure_task_example.R                   Figure 1: task diagram
│   ├── figure_m3_activations.R                 Figure 2: activation decomposition
│   └── 00_download_osf.R                       Download fitted models and figures from OSF
│
├── data/
│   ├── Oberauer_2019_SimpleSpan_Exp1.dat       Tutorial 1: trial-level data
│   ├── Oberauer_2019_SimpleSpan_Exp2.dat       Tutorial 1: trial-level data
│   ├── Oberauer_2019_SimpleSpan_agg.csv        Tutorial 1: aggregated
│   ├── Li_2026_ComplexSpan_Exp1.csv            Tutorials 2 and 3: trial-level data
│   └── Li_2026_ComplexSpan_Exp1_agg.csv        Tutorials 2 and 3: aggregated
│
└── functions/
    └── clean_plot.R                            Plot theme and colour palette
```

`output/` (fitted models, `.rds`) and `figures/` are not tracked because of their size. They are archived on [OSF](https://osf.io/yb7wm/). The long-run fits with 80,000 posterior draws (about 3 GB) are not archived; the scripts refit them if needed.

## Requirements

- R ≥ 4.1 and [CmdStan](https://mc-stan.org/cmdstanr/) (analyses used CmdStan 2.38 via cmdstanr 0.9.0)
- bmm ≥ 1.3.2 (CRAN). Later bmm versions may change default priors, so refits with them can differ slightly.
- brms 2.23.0

```r
install.packages("pacman")
pacman::p_load(here, bmm, brms, tidyverse, tidybayes, patchwork, loo, bridgesampling, osfr)
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

## Reproducing the analyses

Open `tutorial-m3-bmm.Rproj` so that `here()` finds the project root. There are two ways to reproduce the results.

**Use the archived fits.** Run `scripts/00_download_osf.R`. It downloads `output/`, `figures/`, and `data/` from OSF. The scripts pass `file =` to `bmm()`, so they load a cached fit instead of refitting whenever the file exists. You can then run any tutorial script or render the manuscript.

**Refit everything.** Run the scripts in this order. Some read the fits of earlier scripts.

1. `prepare_Li2026_data.R` (optional; the aggregated data are included)
2. `tutorial1_simple_span.R`, which includes the long-run fits for bridge sampling
3. `tutorial2_complex_span.R`
4. `tutorial3_custom_filtering.R` (reads `output/fit_m3_cs.rds` from Tutorial 2)
5. `tutorial4_parameter_recovery.R`
6. `supplement1_loo_choice_rule.R` (reads the Tutorial 1 long-run fits)
7. `supplement1_savage_dickey_stability.R` (reads the default and long-run fits of Tutorials 1 to 3)
8. `supplement2_parameter_recovery_updating_simple.R` and `supplement2_parameter_recovery_updating.R`
9. `figure_task_example.R` and `figure_m3_activations.R`

The long-run fits and the Supplement 2 design grid take several hours.

To render the manuscript, run `quarto render manuscript/tutorial-m3-bmm.qmd` (the supplements work the same way). The PDF uses the vector figures in `figures/`, and the Word version uses their PNG copies.

## Data sources

- **Oberauer (2019)**: Oberauer, K. (2019). Working memory capacity limits memory for bindings. *Journal of Cognition, 2*(1), 40. <https://doi.org/10.5334/joc.86>. Original data: <https://osf.io/qy5sd/>
- **Li et al. (2026)**: Li, C., Frischkorn, G. T., & Oberauer, K. (2026). Can we process information without encoding it into working memory? *Journal of Experimental Psychology: Learning, Memory, and Cognition*. <https://doi.org/10.1037/xlm0001585>. Original data: <https://osf.io/wpcx5/>

Please cite the original articles when reusing these data.

## License

The code (`scripts/`, `functions/`) is released under the MIT License (see `LICENSE`). The manuscript text and the data from Li et al. (2026) are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The Oberauer (2019) data are redistributed from their original OSF repository; please refer to it for terms of reuse.

## References

- Frischkorn, G. T., & Popov, V. (2025). A tutorial for estimating Bayesian hierarchical mixture models for visual working memory tasks: Introducing the Bayesian Measurement Modeling (bmm) package for R. *Behavior Research Methods, 57*(5), 144. <https://doi.org/10.3758/s13428-025-02643-0>
- Oberauer, K., & Lewandowsky, S. (2019). Simple measurement models for complex working-memory tasks. *Psychological Review, 126*(6), 880–932. <https://doi.org/10.1037/rev0000159>
- bmm package: <https://venpopov.github.io/bmm/>
