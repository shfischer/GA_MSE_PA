data-limited catch rule optimisation with MSE and Genetic Algorithms
================

Using machine learning (genetic algorithms) to optimise the data-limited WKMSYCat34 catch rule 3.2.1
====================================================================================================

Introduction
------------

This repository contains the code for optimising WKMSYCat34 catch rule 3.2.1 with genetic algorithms. The simulation is based on the Fisheries Library in R ([FLR](http://www.flr-project.org/)) and the Assessment for All (a4a) standard MSE framework ([`FLR/mse`](github.com/FLR/mse)) developed during the Workshop on development of MSE algorithms with R/FLR/a4a ([Jardim et al., 2017](https://ec.europa.eu/jrc/en/publication/assessment-all-initiativea4a-workshop-development-mse-algorithms-rflra4a)).

Repository structure
--------------------

The root folder contains the following R scripts:

-   `OM.R`: This script creates the operating models,
-   `funs.R` contains functions and methods used for the creation of the operating models and for running the MSE,
-   `GA_funs.R` contains the function used in the optimisation procedure,
-   `run.R` is an R script for running MSE scenarios and is called from a job submission script
-   `run*.qsub` are job submission scripts which are used on a high performance computing cluster and call `run.R`
-   `run_analyse.R` is for analysing the results

R, R packages and version info
------------------------------

The MSE simulation was run on a high performance computing cluster:

``` r
sessionInfo()
R version 3.5.0 (2018-04-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS release 6.7 (Final)
```

The framework uses FLR and the required FLR packages can be installed with `remotes`:

``` r
remotes::install_github(repo = "flr/FLCore")
remotes::install_github(repo = "flr/FLash")
remotes::install_github(repo = "flr/FLife", ref = "ff95b029e7af5fe72332321a84f2a5f32bcf1850")
```

Furthermore, a fork of the `mse/flr` package is required:

``` r
remotes::install_github(repo = "shfischer/mse", ref = "mseDL")
```

And a modified version of the `GA` package for genetic algorithms which also runs on HPCs:

``` r
remotes::install_github(repo = "shfischer/GA")
```

Furthermore, some more R packages available from CRAN are required:

``` r
install.packages(c("foreach", "DoParallel", "dplyr", "tidyr", "data.table")) 
```
