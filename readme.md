Using a genetic algorithm to optimise a data-limited catch rule
================

## Introduction

This repository contains the code for the publication:

> Fischer, S. H., De Oliveira, J. A. A., Mumford, J. D., and Kell, L. T.
> (2021). Using a genetic algorithm to optimize a data-limited catch
> rule. ICES Journal of Marine Science.
> <https://dx.doi.org/10.1093/icesjms/fsab018>.

The status of the repository used for the publication is saved in the
release [code for ICESJMS
publication](https://github.com/shfischer/GA_MSE/releases/tag/v1.0) and
has a DOI:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4475589.svg)](https://doi.org/10.5281/zenodo.4475589)

The simulation is based on the Fisheries Library in R
([FLR](http://www.flr-project.org/)) and the Assessment for All (a4a)
standard MSE framework ([`FLR/mse`](github.com/FLR/mse)) developed
during the Workshop on development of MSE algorithms with R/FLR/a4a
([Jardim et al.,
2017](https://ec.europa.eu/jrc/en/publication/assessment-all-initiativea4a-workshop-development-mse-algorithms-rflra4a)).

The operating models provided as an input are those from the repository
[shfischer/wklifeVII](https://github.com/shfischer/wklifeVII) as
described in:

> Fischer, S. H., De Oliveira, J. A. A., and Laurence T. Kell (2020).
> Linking the performance of a data-limited empirical catch rule to
> life-history traits. ICES Journal of Marine Science, 77: 1914-1926.
> <https://doi.org/10.1093/icesjms/fsaa054>.

The [PA branch](https://github.com/shfischer/GA_MSE/tree/PA) (also
available from
[shfischer/GA\_MSE\_PA](https://github.com/shfischer/GA_MSE_PA)) of this
repository includes the optimisation with specific risk limits for the
ICES precautionary approach (PA) and contains the code for the
publication:

> Fischer, S. H., De Oliveira, J. A. A., Mumford, J. D., and Kell, L. T.
> (in press). Application of explicit precautionary principles in
> data-limited fisheries management. ICES Journal of Marine Science.

## Repository structure

The root folder contains the following R scripts:

-   `OM.R`: This script creates the operating models,
-   `funs.R` contains functions and methods used for the creation of the
    operating models and for running the MSE,
-   `funs_GA.R` contains the function used in the optimisation
    procedure,
-   `run_ms.R` is an R script for running MSE projections and is called
    from a job submission script
-   `run*.pbs` are job submission scripts which are used on a high
    performance computing cluster and call `run_ms.R`
-   `analysis.R` is for analysing the results

The following input files are provided:

-   `input/stocks.csv` contains the stock definitions and life-history
    parameters
-   `input/brps.rds` contains the FLBRP objects which are the basis for
    the OMs

The following outputs summarising the results from running the
optimisation are provided:

-   `output/pol_obj_fun_explorations_stats.csv` exploration of fitness
    functions for pollack
-   `output/pol_interval_MSY_stats.csv` impact of fixing the catch
    advice interval for pollack
-   `output/all_stocks_MSY_stats.csv` optimisation results for all 29
    simulated stocks
-   `output/groups_MSY_stats.csv` optimisation results for stock groups

## R, R packages and version info

The MSE simulations were run on a high performance computing cluster:

``` r
> sessionInfo()
R version 3.6.1 (2019-07-05)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /rds/general/user/shf4318/home/anaconda3/envs/R_2020/lib/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] doMPI_0.2.2         Rmpi_0.6-9          doRNG_1.8.2
 [4] rngtools_1.5        doParallel_1.0.15   GA_3.2.1
 [7] foreach_1.4.8       mse_2.0.3           FLBRP_2.5.4
[10] data.table_1.12.2   ggplotFL_2.6.7.9001 ggplot2_3.1.1
[13] FLash_2.5.11        FLCore_2.6.14.9004  iterators_1.0.12
[16] lattice_0.20-40

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.5       pillar_1.4.6     compiler_3.6.1   plyr_1.8.4
 [5] tools_3.6.1      digest_0.6.18    lifecycle_0.2.0  tibble_2.1.1
 [9] gtable_0.3.0     pkgconfig_2.0.2  rlang_0.4.5      Matrix_1.2-18
[13] cli_2.0.2        gridExtra_2.3    withr_2.3.0      dplyr_0.8.0.1
[17] stats4_3.6.1     grid_3.6.1       tidyselect_0.2.5 glue_1.3.2
[21] R6_2.4.0         fansi_0.4.1      purrr_0.3.3      magrittr_1.5
[25] scales_1.0.0     codetools_0.2-16 ellipsis_0.3.0   MASS_7.3-51.5
[29] assertthat_0.2.1 colorspace_1.4-1 lazyeval_0.2.2   munsell_0.5.0
[33] crayon_1.3.4
```

The framework is based on the Fisheries Library in R (FLR) framework.
The exact versions of the packages as used here can be installed with
`remotes`:

``` r
remotes::install_github(repo = "flr/FLCore", ref = "3d694903b9e6717b86c3e8486fc14ebf92908786")
remotes::install_github(repo = "shfischer/FLash", ref = "d1fb86fa081aaa5b6980d74b07d9adb44ad19a7f", INSTALL_opts = "--no-multiarch") # silenced version of FLash
# INSTALL_opts = "--no-multiarch" to avoid issues in Windows
remotes::install_github(repo = "flr/FLBRP", ref = "3a4d6390abc56870575fbaba3637091036468217", INSTALL_opts = "--no-multiarch")
```

Furthermore, a data-limited fork of the `flr/mse` package is required:

``` r
remotes::install_github(repo = "shfischer/mse", ref = "mseDL2.0", INSTALL_opts = "--no-multiarch")
```

And a modified version of the `GA` package for genetic algorithms which
also runs on HPCs and supports MPI parallelisation:

``` r
remotes::install_github(repo = "shfischer/GA")
```

Furthermore, some more R packages available from CRAN are required:

``` r
install.packages(c("foreach", "DoParallel", "doRNG", "dplyr", "tidyr", "ggplot2", "scales", "cowplot", "Cairo", "scales")) 
```

For using MPI parallelisation, an MPI backend such as OpenMPI and the R
packages `Rmpi` and `doMPI` are required.
