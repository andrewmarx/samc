## **samc** <a href="reference/figures/3d-stack.png"><img align="right" width=35% src="man/figures/3d-stack-small.png" style="padding-left: 10px"></a>

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/samc)](https://cran.r-project.org/package=samc)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/samc)](https://cran.r-project.org/package=samc)


#### **Introduction**

This is an R package that implements functions for working with absorbing Markov chains using theorems described in the book "Finite Markov Chains" by Kemeny and Snell. The design of this package is based on the implementation of these theorems in the framework described by Fletcher et al. (2019), which applies them to spatial ecology. This framework incorporates both resistance and absorption (or mortality) using spatial absorbing Markov chains to provide several short- and long-term predictions for metrics related to connectivity in landscapes. These metrics are listed in Table 1 of Fletcher et al. (2019), as well as the [Overview](articles/overview.html) vignette. Despite the ecological context of the package, these functions can be used in any application of absorbing Markov chains.

As of version 3.1.0, this package also supports use of the convolution algorithm described by Hughes et al. (2023). This algorithm uses a moving window approach to calculate many of the same metrics iteratively in a fast and memory efficient manner.


#### **Installation**

It is recommended that users install the samc package via CRAN, where it will be regularly kept up to date.

```R
install.packages("samc")
```


#### **Citation**

When using the package, please cite it using:

Marx, A.J., Wang, C., Sefair, J.A., Acevedo, M.A. and Fletcher, R.J., Jr. (2020), samc: an R package for connectivity modeling with spatial absorbing Markov chains. Ecography, 43: 518-527. [doi:10.1111/ecog.04891](https://doi.org/10.1111/ecog.04891)


#### **Publications**

This is a list of publications where the samc package was used in some way:

- [**Predicting dispersal and conflict risk for wolf recolonization in Colorado**](http://dx.doi.org/10.1111/1365-2664.14504) *Journal of Applied Ecology (Sep 2023)*
- [**A framework for linking dispersal biology to connectivity across landscapes**](http://dx.doi.org/10.1007/s10980-023-01741-8) *Landscape Ecology (Jul 2023)*
- [**Comparison and parallel implementation of alternative moving-window metrics of the connectivity of protected areas across large landscapes**](http://dx.doi.org/10.1007/s10980-023-01619-9) *Landscape Ecology (Mar 2023)*
- [**Mapping the connectivity–conflict interface to inform conservation**](http://dx.doi.org/10.1073/pnas.2211482119) *PNAS (Dec 2022)*
- [**Landscape connectivity for an endangered carnivore: habitat conservation and road mitigation for ocelots in the US**](http://dx.doi.org/10.1007/s10980-022-01569-8) *Landscape Ecology (Dec 2022)*
- [**Extending isolation by resistance to predict genetic connectivity**](http://doi.org/10.1111/2041-210X.13975) *Methods in Ecology and Evolution (Sep 2022)*
- [**Defining and quantifying effective connectivity of landscapes for species' movements**](http://dx.doi.org/10.1111/ecog.05351) *Ecography (Mar 2021)*


#### **Major Changes**

Version 4 of the package made one breaking change:

- The `visitation_net(samc, origin, dest)` function has been updated to behave consistently like other metrics. Instead of a returning a vector, it now returns a single value. This is accompanied with a change to the mathematical implementation the correctly generalizes to models with multiple absorbing points.

Version 3 of the package made some minor breaking changes:

- The `samc()` function no longer supports `TransitionLayer` inputs. This only had a niche use case, but before v3 wasn't an issue to include because of other dependencies on gdistance. With v3, this became the only dependency left for gdistance, so it was removed to avoid potential future issues should gdistance ever get removed from CRAN (which nearly happened in 2022).
- With the addition of terra support, the `map()` function was updated so that its output matches the input type to `samc()`. Previously, matrix inputs were matched to RasterLayers, but now they are mapped back to matrices.
- The `sym` option for creating the samc object is currently ignored.
- The `tr_args` and `occ` parameters were renamed to `model` and `init`, respectively.
- `cond_passage()` and `visitation()` had an `init` argument inserted to match the usage of other metrics. These arguments are not fully implemented as of v3.0.0 but may be implemented in the future.
- Cells are no longer automatically named by `samc()` when creating the transition matrix from maps. Generating unique character names for each cell ended up being a significant waste of memory as the inputs to `samc()` got larger. It's unlikely these names were ever used in practice, since using the numeric results from `locate()` is more convenient.

Version 2 of the package officially removed support for various deprecated parameters in the `samc()` function. Deprecation warnings were provided starting in v1.4.0 of the package, along with message details and a backward compatible implementation of the expected changes. Removing this backward compatibility is a breaking change that will require some old code to be updated in order to run on the latest version of the package. The changes needed are straightforward and mostly entail some reorganization of the input parameters for the `samc()` function. Some of the old functionality, primarily overriding memory safety limits, has been moved to the `samc-class` itself and is no longer tied to the object creation. Redesigning the `samc()` function and removing backward compatibility makes maintaining the package and adding new features later a substantially improved process; hopefully with only minor inconvenience to users.


#### **Support**

Please note that this section is for package-specific queries. If you have questions or comments about the related published articles, contact the authors.

###### *Have an error and can't fix it?*

If you have an error, make sure you are using the newest version of the package. Then, check the site's [Troubleshooting](articles/troubleshooting.html) page, which is periodically updated with common errors that users encounter.

###### *Not sure how to do something? Found a bug? Have a feature request? Still can't solve your error? Want to show off your project?*

We have a Github discussions page for anything and everything related to the package here: [Github Discussions](https://github.com/andrewmarx/samc/discussions)


#### **License**

By default, package code is licensed under AGPLv3. Some code files derived from 3rd party sources may be licensed under GPLv3 and will have a comment to indicate this is the case.

Website materials are licensed under CC BY-NC-SA 4.0.


#### **References**

Fletcher, R.J., Jr., Sefair, J.A., Wang, C., Poli, C.L., Smith, T.A.H., Bruna, E.M., Holt, R.D., Barfield, M., Marx, A.J. and Acevedo, M.A. (2019), Towards a unified framework for connectivity that disentangles movement and mortality in space and time. Ecol Lett, 22: 1680-1689. https://doi.org/10.1111/ele.13333

Hughes, J., Lucet, V., Barrett, G. et al. Comparison and parallel implementation of alternative moving-window metrics of the connectivity of protected areas across large landscapes. Landsc Ecol (2023). https://doi.org/10.1007/s10980-023-01619-9
