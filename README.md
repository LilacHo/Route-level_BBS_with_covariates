# Route-level BBS with covariates

Estimate route-level population trends from the North American Breeding Bird
Survey (BBS) while accounting for land-cover change, then compare those trends
against Species Distribution Model (SDM) projected range-change categories.

For each BBS route the model separates the full population trend from a
residual trend net of land-cover change (`CH_no_habitat`). Route-level residual
trends are then tested across SDM change categories (RCP 4.5 and RCP 8.5) and
visualized. The modelling approach builds on 
[Smith et al. 2024](https://ace-eco.org/vol19/iss2/art23/).

## The `land_cover` parameter

Scripts 1–3 share a single `land_cover` setting (e.g. `"grasslands"`) near the
top of each file. It must match a value in the `Group` column of
`data/spp_names_codes_group.csv` and drives all input and output names:
the covariate file (`data/<land_cover>.csv`), the SDM raster folders
(`data/rcp45_<land_cover>/`, `data/rcp85_<land_cover>/`), and the
`<land_cover>_`-prefixed output files. Change the one value in each script to
run a different species group.

## Pipeline

Run the scripts in `code/` in order. Each one consumes the output of the
previous step.

| Step | Script | Purpose |
|------|--------|---------|
| 0 | `0_prepare_aou.R` | Add AOU codes / BBS names to the species list |
| 1 | `1_CH_no_habitat_routes.R` | Fit the trend model; compute per-route `CH`, `CH_no_habitat`, `CH_dif` |
| 2 | `2_add_SDM.R` | Add SDM range-change categories (`rcp45`, `rcp85`) to each route |
| 3 | `3_statistical_analysis_and_visualization.R` | Test and plot residual trends by SDM category |

`functions/neighbours_define_voronoi.R` builds the spatial adjacency used by
the model; `models/slope_habitat_route_NB.stan` is the fitted Stan model.

### Step 0 — `0_prepare_aou.R`
Adds AOU numeric codes (via `wildlifeR`) to the source species list and matches
each species to the `bbsBayes2` species list, recording the exact `bbs_english`
name used to query BBS.

- **Input:** `data/spp_names_codes_group.csv` (available from [Bateman et al. 2020 data](https://adaptwest.databasin.org/pages/audubon-survival-by-degrees/))
- **Output:** `data/spp_names_codes_group_aou.csv`

### Step 1 — `1_CH_no_habitat_routes.R`
For every species in the target group, pulls BBS counts (2010–2024 by default),
keeps routes that fall within the strata and have
covariate data, builds spatial neighbours, and fits the negative-binomial
trend model with the land-cover covariate. Converts the full and residual
posterior slopes into cumulative percent change per route. Fits and Stan data
are cached, and existing per-species CSVs are skipped, so the script is
resumable.

The model and data-prep code for this step are adapted from
[AdamCSmithCWS/Route-level_BBS_trends](https://github.com/AdamCSmithCWS/Route-level_BBS_trends)
and [AdamCSmithCWS/Jefferys_etal](https://github.com/AdamCSmithCWS/Jefferys_etal).

- **Input:** `data/spp_names_codes_group_aou.csv`, `data/<land_cover>.csv`,
  `models/slope_habitat_route_NB.stan`, BBS data fetched via `bbsBayes2`
- **Output:** `output/species_routes/<land_cover>_<species>_route_CH.csv`
  (columns: `species`, `species_code`, `land_cover`, `route`, `routeF`,
  `latitude`, `longitude`, `CH`, `CH_no_habitat`, `CH_dif`); cached fits in
  `output/` and `data/stan_data/`, route maps in `data/maps/`

### Step 2 — `2_add_SDM.R`
For each route CSV, extracts the SDM classified-change raster value
([Bateman et al. 2020](https://conbio.onlinelibrary.wiley.com/doi/10.1111/csp2.242);
[data access](https://adaptwest.databasin.org/pages/audubon-survival-by-degrees/))
at the route location for both climate scenarios and adds them as `rcp45` and
`rcp85` columns. Missing rasters produce a warning and `NA` values.

- **Input:** `output/species_routes/*_route_CH.csv`,
  `data/rcp45_<land_cover>/<CODE>/...tif`,
  `data/rcp85_<land_cover>/<CODE>/...tif`,
  `data/spp_names_codes_group_aou.csv`
- **Output:** `output/species_routes_sdm/<land_cover>_<species>_route_CH_sdm.csv`

Raster value legend: `0` never suitable, `1` extirpation,
`2` worsening, `3` slightly worsening, `4` neutral, `5` slightly improving,
`6` improving, `7` colonization.

### Step 3 — `3_statistical_analysis_and_visualization.R`
Tests whether residual trends (`CH_no_habitat`) differ across SDM categories —
using Kruskal–Wallis plus pairwise Wilcoxon (BH-adjusted) tests with compact
letter displays — and draws faceted RCP 4.5 / 8.5 violin plots. Runs both the
full 8-category scheme (0–7) and a grouped scheme (Contraction = 1,2,3 /
Stable = 4 / Expansion = 5,6,7; category 0 excluded), each for all species
combined and per species. The analysis is run twice: first on the selected
`land_cover` only, then on **all** land covers pooled together (every
`*_route_CH_sdm.csv` in the folder, regardless of prefix).

- **Input:** `output/species_routes_sdm/<land_cover>_*_route_CH_sdm.csv` (for
  the selected group), plus all other `*_route_CH_sdm.csv` for the pooled
  all-land-cover analysis
- **Output:**
  - stats reports in `output/species_routes_sdm_stats/`
  - all-species violin plots in `output/species_routes_sdm_plot/`
  - per-species violin plots in `output/species_routes_sdm_plot/per_species/`
  - selected-group outputs are prefixed with `<land_cover>_`; the pooled
    outputs are prefixed with `all_lc_`

## Data

- `data/spp_names_codes_group_aou.csv` — species lookup: `Common.Name`, `Code`
  (4-letter alpha), `Group`, AOU `spp.num`, `bbs_english`. Produced by Step 0.
- `data/<land_cover>.csv` — route-level land-cover covariate by year (e.g.
  `grasslands.csv`, `developed.csv`); route key is `StateNum-Route`.
- `data/rcp45_<land_cover>/<CODE>/...classifiedchange.tif` and
  `data/rcp85_<land_cover>/<CODE>/...` — SDM classified-change rasters per
  species, for the RCP 4.5 and RCP 8.5 scenarios.

The route-level land-cover covariate (`data/<land_cover>.csv`) can be generated
by following [LilacHo/BBS_landcover](https://github.com/LilacHo/BBS_landcover).

## Requirements

- R (≥ 4.1 recommended), plus the matching
  [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (Windows) or
  Xcode command-line tools (macOS) — a C++ toolchain is needed to compile
  CmdStan in Step 1. Install the Rtools version that pairs with your R version.
- [`cmdstanr`](https://mc-stan.org/cmdstanr/) and a working CmdStan install
  (required for Step 1):

  ```r
  install.packages("cmdstanr",
                   repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
  cmdstanr::install_cmdstan()
  ```

- `bbsBayes2`:

  ```r
  install.packages("bbsBayes2",
                   repos = c(bbsbayes = 'https://bbsbayes.r-universe.dev',
                             CRAN = 'https://cloud.r-project.org'))
  ```

  After installing, run `bbsBayes2::fetch_bbs_data()` the first time you use it
  to download the BBS data locally. This only needs to be done once — subsequent
  runs reuse the downloaded data.

- Remaining R packages: `here`, `tidyverse`, `posterior`, `sf`, `terra`,
  `multcompView`, `wildlifeR`:

  ```r
  install.packages(c("here", "tidyverse", "posterior", "sf", "terra",
                     "multcompView", "wildlifeR"))
  ```




