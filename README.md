# samc <a href="reference/figures/3d-stack.png"><img align="right" width=35% src="man/figures/3d-stack-small.png" style="padding-left: 10px"></a>

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/samc)](https://cran.r-project.org/package=samc)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/samc)](https://cran.r-project.org/package=samc)
[![Gitter](https://badges.gitter.im/samc-package/community.svg)](https://gitter.im/samc-package/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)


## Introduction

This is an R package that implements functions for working with absorbing Markov chains using theorems described in the book "Finite Markov Chains" by Kemeny and Snell. The design of this package is based on the implementation of these theorems in the framework described by Fletcher et al. in [*Toward a unified framework for connectivity that disentangles movement and mortality in space and time*](https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.13333) (Ecology Letters, 2019; DOI: [10.1111/ele.13333](https://doi.org/10.1111/ele.13333)), which applies them to spatial ecology. Despite the ecological context of the package, these functions can be used in any application of absorbing Markov chains.

For more information, and recommended installation instructions, most users should visit the main home page at https://andrewmarx.github.io/samc. The remainder of this document is dedicated to installing the package from source (not recommended).


## Installing From Github

#### Pre v1.0.0

Specific versions of the package can be installed by using the `ref` parameter. In the case of Fletcher et al. (2019), version 0.1.0 of the package was used. To run the code example in the SI of that paper, you would need to install that specific version of samc using:

```R
devtools::install_github("andrewmarx/samc", ref = "0.1.0")
```

#### v1.0.0 and Later

Version 1.0.0 and newer requires C++ development tools in order to install from source. The steps required to install the appropriate development tools varies by operating system and is beyond the scope of this document. Users will have to locate and follow appropriate external documentation to setup the devtools if they wish to install the newest version of the samc package from source.

If the devtools are installed and setup correctly, then the latest version of the package can be installed directly from github using the following command:

```R
devtools::install_github("andrewmarx/samc")
```
