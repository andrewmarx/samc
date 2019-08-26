# samc 1.0.3

- First CRAN submission


# samc 1.0.0

- Complete package rewrite (code dependent on v0.1.0 will not work)
- samc-class for managing SAMC data
- Utility functions for creating samc-class objects, checking inputs, and mapping data
- Heavily optimized analytical functions for all metrics described in Fletcher et al. (2019, DOI: [10.1111/ele.13333](https://doi.org/10.1111/ele.13333))
  - Utilizing sparse matrices
  - Eigen C++ implementation via Rcpp and RcppEigen
- Updated example data
- Extensive documentation
- Several tutorials


# samc 0.1.0

- Created crude functions for calculating metrics in Fletcher et al. (2019, DOI: [10.1111/ele.13333](https://doi.org/10.1111/ele.13333))
- Included example data
