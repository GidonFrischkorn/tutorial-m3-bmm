# Codebook for `data/`

Variable keys for the five data files used in the tutorials. Column meanings for the
trial-level files follow the original analysis scripts and data repositories linked under
each source. Where this repository does not document a variable's meaning, the entry says so
and points to the original repository.

## Oberauer (2019): simple span (Tutorial 1)

Source: Oberauer, K. (2019). Working memory capacity limits memory for bindings. *Journal of
Cognition, 2*(1), 40. <https://doi.org/10.5334/joc.86>. Original data: <https://osf.io/qy5sd/>.

### `Oberauer_2019_SimpleSpan_Exp1.dat`, `Oberauer_2019_SimpleSpan_Exp2.dat`

Trial-level data. Space-delimited, **no header row**. Exp. 1 is the open-set experiment,
Exp. 2 the closed-set experiment (between subjects). Participant IDs overlap across the two
files. Column names below are those of the original JAGS analysis script
(<https://osf.io/qy5sd/files/geb6v>), in file order:

| # | Column | Description |
| --- | --- | --- |
| 1 | `id` | Participant ID (unique within an experiment only) |
| 2 | `session` | Session (1–3) |
| 3 | `block` | Block within session (1–8) |
| 4 | `trial` | Trial within block (1–28) |
| 5 | `setsize` | Number of studied items (2, 4, 6, 8) |
| 6 | `rsizeList` | Number of studied items in the response set, including the correct item (0 = recall trial) |
| 7 | `rsizeNPL` | Number of not-presented lures in the response set (0 = recall trial) |
| 8 | `tested` | Probed serial position (1–8) |
| 9 | `response` | Code of the selected response item (stimulus code; see original repository) |
| 10 | `rcat` | Response category: 1 = correct, 2 = other list item, 3 = not-presented lure |
| 11 | `rt` | Response time in seconds; −1 on recall trials (no response time recorded) |

The tutorial uses recognition trials only (`rsizeList + rsizeNPL > 0`); see
`scripts/tutorial1_simple_span.R`, Section 0.2.

### `Oberauer_2019_SimpleSpan_agg.csv`

Aggregated by `scripts/tutorial1_simple_span.R` (Sections 0.1–0.7). One row per participant ×
recognition condition (960 rows: 40 participants × 24 conditions).

| Column | Description |
| --- | --- |
| `id` | Participant ID, made unique as `<id>_<exp>` |
| `exp` | Experiment: `openset` (Exp. 1) or `closedset` (Exp. 2) |
| `setsize` | Set size (2, 4, 6, 8) |
| `ss_fac` | Set size as a factor |
| `ss_lin` | Set size centered on its mean (linear predictor) |
| `corr`, `other`, `npl` | Response frequencies: correct, other list item, not-presented lure |
| `n_corr`, `n_other`, `n_npl` | Number of response options in each category (`n_corr` = 1; `n_other` = `rsizeList` − 1; `n_npl` = `rsizeNPL`) |

## Li, Frischkorn, & Oberauer (2026): complex span (Tutorials 2 and 3)

Source: Li, C., Frischkorn, G. T., & Oberauer, K. (2026). Can we process information without
encoding it into working memory? *Journal of Experimental Psychology: Learning, Memory, and
Cognition*. <https://doi.org/10.1037/xlm0001585>. Original data: <https://osf.io/wpcx5/>.

### `Li_2026_ComplexSpan_Exp1.csv`

Trial-level data of Experiment 1, one row per screen (18,000 rows: 50 participants × 60 trials
× 6 screens). Each trial shows three image pairs with a size judgment (memory screens); one
image per pair is cued as the memory target. At retrieval, participants select each target
from 12 images.

| Column | Description |
| --- | --- |
| `participant` | Participant ID |
| `blockID` | Block (1–3; one condition per block) |
| `trialID` | Trial within block (1–20) |
| `condition` | `control` (no distractors), `pre` (target cued before the size judgment), `retro` (target cued after the size judgment) |
| `screenID` | `memory` (image pair with size judgment) or `retrieval` (test screen) |
| `leftImage`, `rightImage` | Image IDs shown on a memory screen (`NA` = image failed to load) |
| `cue` | Side of the cued memory target on a memory screen (`left`, `right`) |
| `congruency` | `congruent` / `incongruent` trial label from the original experiment (see original repository) |
| `correctObject` | Image ID of the tested target on a retrieval screen |
| `response` | Memory screen: size-judgment key (`left`, `right`); retrieval screen: selected image ID; `NULL` = no response |
| `acc` | Accuracy of the response on that screen (`TRUE`/`FALSE`) |
| `rt` | Response time in milliseconds |

### `Li_2026_ComplexSpan_Exp1_agg.csv`

Aggregated by `scripts/prepare_Li2026_data.R` after the participant exclusions described there.
One row per participant × condition (135 rows: 45 participants × 3 conditions).

| Column | Description |
| --- | --- |
| `participant` | Participant ID |
| `condition` | `control`, `pre`, `retro` |
| `corr` | Frequency of selecting the tested target |
| `other` | Frequency of selecting another target from the same trial |
| `distc` | Frequency of selecting the distractor paired with the tested target |
| `disto` | Frequency of selecting another distractor from the same trial |
| `npl` | Frequency of selecting a not-presented lure |
| `n_corr`, `n_other`, `n_distc`, `n_disto`, `n_npl` | Number of response options per category among the 12 retrieval images (control: 1, 2, 0, 0, 9; pre/retro: 1, 2, 1, 2, 6) |

## License

The Li et al. (2026) data are released under CC BY-NC 4.0. The Oberauer (2019) data are
redistributed from their original OSF repository; refer to it for terms of reuse. Please
cite the original articles when reusing these data.
