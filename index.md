# samc <a href="reference/figures/3d-stack.png"><img align="right" width=35% src="man/figures/3d-stack-small.png" style="padding-left: 10px"></a>

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/samc)](https://cran.r-project.org/package=samc)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/samc)](https://cran.r-project.org/package=samc)

## Introduction

This is an R package that implements functions for working with the framework described by Fletcher et al. in [*Toward a unified framework for connectivity that disentangles movement and mortality in space and time*](https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.13333) (Ecology Letters, 2019; DOI: [10.1111/ele.13333](https://doi.org/10.1111/ele.13333)). This framework incorporates both resistance and absorption (or mortality) using spatial absorbing Markov chains to provide several short- and long-term predictions for metrics related to connectivity in landscapes. These metrics are listed in Table 1 of Fletcher et al. (2019), as well as the [Overview](articles/overview.html) vignette.


## Installation

It is recommended that users install the samc package via CRAN, where it will be regularly kept up to date.

```R
install.packages("samc")
```


## Citation

Marx, A.J., Wang, C., Sefair, J.A., Acevedo, M.A. and Fletcher, R.J., Jr. (2020), samc: an R package for connectivity modeling with spatial absorbing Markov chains. Ecography, 43: 518-527. [doi:10.1111/ecog.04891](https://doi.org/10.1111/ecog.04891)


## Support

Please note that this section is for package specific queries. If you have questions or comments about the related published articles, contact the authors.

#### Not sure how to do something? Have an error and can't fix it?

If you have an error, make sure you are using the newest version of the package. Then, check the site's [Troubleshooting](articles/troubleshooting.html) page, which is periodically updated with common errors that users encounter. If you still cannot solve your issue, or have another programming related question, then [Stack Overflow](https://stackoverflow.com/) is the place to go. Please respect the rules and etiquette expected there. This includes first checking to see if anyone else has asked similar questions about the samc package. If you can't find an answer, then ask a question. Be sure to include the following information: operating system, R version, samc package version. Also make sure you post relevant code, ideally a self-contained reproducible example that others can run.

#### Found a bug? Have a feature request?

If you're getting an error or incorrect output, and are sure it is a bug, then report it as an issue in the GitHub repo (after making sure that it hasn't already been posted). Make sure to include all possible relevant info (operating system, R version, samc version) and a fully reproducible example. Likewise, if you have a feature you'd like to see, submit it to the issue tracker. Keep in mind that my time is limited, so while I do have plans to add new features, I also have other obligations.
