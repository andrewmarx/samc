
# .onLoad <- function(libname, pkgname) {
#   assign("samc_env", new.env(parent = emptyenv()), parent.env())
# }

.onAttach <- function(libname, pkgname) {
  msg <-paste("Breaking changes in version 3:",
              "\n1) The samc() function no longer supports TransitionLayer inputs.",
              "\n2) The `map()` function will now return matrices if matrices were used to create the samc object.",
              "\n3) `cond_passage()` and `visitation()` had an `occ` argument inserted to match the usage of other metrics.",
              "\n4) Cells are no longer automatically named by `samc()` when creating the transition matrix from maps.",
              "\n",
              "The following parameters were removed from the samc() function in version 2:",
              "\n1) `latlon`: Now handled automatically.",
              "\n2) `resistance`: `data` should be used in its place.",
              "\n3) `p_mat`: The P matrix is now provided via `data`.",
              "\n4) `tr_fun`: See the samc() function documentation for the alternative.",
              "\n5) `directions`: See the samc() function documentation for the alternative.",
              "\n6) `override`: See the samc-class documentation for the alternative.")
  packageStartupMessage(msg)
}
