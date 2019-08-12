# SAMC v0.1.0

## Introduction

This is an R package that implements functions for working with spatial absorbing Markov chains based on the framework described in *Toward a unified framework for connectivity that disentangles movement and mortality in space and time* by Fletcher et al. in Ecology Letters (2019; DOI: [10.1111/ele.13333](https://doi.org/10.1111/ele.13333)). This framework incorporates both resistance and absorption (or mortality) to provide several short- and long-term predictions for metrics related to connectivity in landscapes.

Note: v1.0.0 of this package is currently under development and is expected to be released soon (August, 2019). It includes a complete overhaul of the package, support for extremely large landscapes (1M+ cells), substantial speed and memory improvements, and extensive documentation.

## Installation

### Installing From Github

The most up to date version of the package can be installed directly from github using the following command:
```R
devtools::install_github("andrewmarx/samc")
```

Specific versions of the package can be installed by using the `ref` parameter:
```R
devtools::install_github("andrewmarx/samc", ref = "0.1.0")
```
