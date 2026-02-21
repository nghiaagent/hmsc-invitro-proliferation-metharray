# Genomic Profiling of Drivers of hMSC in vitro Proliferation and Differentiation Potential


# Introduction

This repository contains the source code for DNA methylation analysis in
the thesis project: Genomic Profiling of Drivers of hMSC *in vitro*
Proliferation and Differentiation Potential.

See the thesis here: <https://eprints.qut.edu.au/262996/>

# Prerequisites

Make sure the following are installed on your system:

- Quarto
- R 4.5.1

# Reproduciblity

To reproduce the thesis figures:

- Install all prerequisites software.
- Clone the repository; Positron is recommended for this purpose.
- Request a copy of our dataset via
  <minhnghia.nguyen@connect.qut.edu.au> and place it in the `input`
  folder as instructed.
- Run the following in an R terminal to download all packages:

``` r
renv::restore()
```

- Run the following in an R terminal to perform all analyses:

``` r
source(here::here("R/99_pipeline_firstrun.R"))
```
