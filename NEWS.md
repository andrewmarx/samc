# samc 1.5.0

- Updated `cond_passage()` to return `0` for when *i==j* in the vectors. This fixes an issue associated with shifted indices in `cond_passage(samc, dest)`. It also technically breaks backwards compatibility for when `dest` equals `origin` in `cond_passage(samc, origin, dest)`. Previously, `cond_passage(samc, origin, dest)` would return `NA` when `origin` equaled `dest`, but this decision was arbitrary. The `cond_passage()` documentation explains why.
- Added a new section for worked examples on the website.
- Added a new example illustrating how to use various aspects of the package with a simple perfect maze and interpret the results. See the Maze Part 1 vignette.
- Added multithreading for the `dispersal(samc, origin/occ)` function via the RcppThread package. See the Parallel Computing vignette for details.

# samc 1.4.1

- Added an input check for multiple absorption that throws a more informative error when a list contains anything other than matrices
- Updated the crs check in `samc()`. CRS objects have a hidden field that can vary depending on system and software versions, and previous versions of the check would not account for this. This would to lead to false positives where perfectly compatible rasters were reported as incompatible. The corresponding error message was also fixed to report the correct issue; the code was initially copied and modified from another input check, but the error message wasn't updated in the process.
- Added an initial vignette discussing *Disconnected Data*. The current contents are only slightly modified from an email discussion; they will be rewritten and expanded upon in the future. The *Troubleshooting* vignette has had an error message and a warning message related to the topic added to it.
- Added a Rcpp related error to the *Troubleshooting* vignette.
- Bumped version requirements for R to 3.6.0, Rcpp to 1.0.5, RcppEigen to 0.3.3.9.1, and set C++14 as the standard to use in Makevars.
- Enabled the Github discussions page as a replacement for Gitter: https://github.com/andrewmarx/samc/discussions

# samc 1.4.0

- Due to a ballooning parameter count, the samc() function parameters are being adjusted. The new version is samc(data, absorption, fidelity, tr_args). Code using the previous syntax should continue to work (with one rare edge-case as an exception), but backwards compatibility will be remove in version 1.5.0, so old code should be updated. See the samc() function documentation and website tutorials for full details and examples. Package startup output has been added to detail the changes as well.
  - Updated long/lat handling in samc() to use projection info built into the raster. Deprecated latlon parameter (no longer needed). Added warning for when rasters have non-square cells and are missing projection information.
  - The data parameter should be used to pass in the data related to transition probabilities (essentially replaces the resistance and p_mat parameters)
  - The tr_fun and directions parameters have been deprecated. This information is now passed as list to the tr_args.
  - Deprecated override parameter in the samc() function. See samc-class documentation for details on how to set this.
- Added support for specifying if transition functions are symmetrical or not through the tr_args parameter list.
- Added the ability to directly input a custom TransitionLayer to the samc() function. This allows more flexibility than RasterLayer/matrix maps, but is a little safer than directly inputting a P matrix. See samc() documentation and *Overview* vignette for full details.
- Added the ability to use the $ operator for accessing and modifying components of samc-class objects. See samc-class documentation for details.
- Updated check() so that multiple rasters can be inputted in the first argument as a RasterStack. This eliminates the need to manually run check() for multiple pairs of rasters.
- Added initial support for caching intermediate results of some calculations. This currently only benefits dispersal(samc, occ), which now caches the diag(F) calculation. This means that while the first run of this method will still be slow, subsequent runs will be substantially faster. With this feature, dispersal(samc, origin) has been enabled, and will share the same cached information with dispersal(samc, occ). Future versions will expand the cache options to additional metrics.
- Added support for multiple absorption. The `absorption` parameter in samc() is treated as the total absorption (consistent with previous behavior). After creation of the samc-class object, additional absorbing states can be attached to the samc-class object. See the samc-class documentation and the new *Multiple Absorption* tutorial for more details. 
- Added a new absorption() metric. This metric is closely related to the mortality() metric. The absorption() metric can be used to determine the overall probability that a particular absorbing state will be reached (the mortality() metric calculates it for individual transient states rather than overall).
- Fix missing value short-term dispersal
- Overhauled the *Overview* vignette, including adding more details about the construction of the P matrix.
- Performance vignette update
- Updated documentation for various analytical functions, including more formal/consistent terminology.
- Vector outputs from metric should now all have named cells. These names correspond to the row/column names of the P matrix.

# samc 1.3.0

- Fixed an issue with the check() function when data contains NA's.
- Fixed an issue with the raster returned from locate(samc) having 0 for NA cells.
- Improved error checking and messaging for the check() and locate() functions.
- Named rows and columns for the P matrix is now supported. Previously, naming the rows and columns would cause some checks to fail. If names are not manually assigned, the names are simply the row/column numbers converted to character strings.
- Analytical functions updated to support named inputs for the origin and dest location parameters
- When both the origin and dest parameter is used in a function, the inputs can be paired vectors.
- Added the pairwise() utility function
- Created a new *Locations* tutorial vignette for new location input options.

# samc 1.2.1

- Fixed a regression in v1.2.0 where the samc() function would not work correctly unless matrix/raster layers contained at least one NA cell
- Revamped the automated test suite with more test scenarios to better catch issues before release
- Added checks during samc-class creation to prevent potential issues with discontinuous/clumped input data. Currently, this type of data will not work with the cond_passage() function, but will in a future release.
- Reworked some of the vignettes to produce cleaner pages and remove suggested dependencies (e.g. gifski, gganimate, ggplot2) from the package so that users aren't bugged about installing them if they don't need them.


# samc 1.2.0

- Added the ability to create samc-class objects from a custom P matrix using p_mat parameter in samc(). See the samc() documentation for details
- Added the cond_passage() function, which calculates conditional mean first passage times
- Added the locate() function, which functions similarly to the cellFromXY() function in the raster package. It's used to get cell numbers from xy coords, but unlike cellFromXY(), it properly accounts for how cells are numbered when the P matrix is constructed.
- Adjusted the absorption inputs to support values of 0 (i.e., no absorption). Currently, at least one cell must have a non-zero value
- Fixed an issue where raster/matrix inputs containing isolated cells (individual cells neighbored by only NA values) would lead to malformed P matrices.


# samc 1.1.0

- Added support for the use vectors of time steps in most short-term metrics. It is more computationally efficient and ergonomic to do this rather than calculating short-term metrics for each individual time step. Some key points:
  - When vector inputs are used for time steps, the result is contained in a list
  - The names of the entries in the list are character versions of the corresponding time step values
  - Time step vectors must consist of ordered positive integers with no duplicate values
  - Time step vector inputs have not been added for short-term metrics that return dense matrices
- Updated the map() function to support list inputs generated by the short-term metrics. The result is a list of RasterLayers
- Updated the *Temporal Analysis* and *Animations* vignettes to incorporate time step vectors
- Created a [Gitter community](https://gitter.im/samc-package/community) for package support. Gitter badges on the README and home pages can be used to access it.
- Updated the package citation info to refer to Marx et al. (2020, DOI: [10.1111/ecog.04891](https://doi.org/10.1111/ecog.04891))


# samc 1.0.4

- Add conditional usage of suggested packages in vignettes
- Minor updates for package info


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
