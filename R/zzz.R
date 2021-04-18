
# .onLoad <- function(libname, pkgname) {
#   assign("samc_env", new.env(parent = emptyenv()), parent.env())
# }

.onAttach <- function(libname, pkgname) {
  msg <-paste("Important changes have been made to the samc() function starting in version 1.4.0:",
              "\n1) `latlon` is no longer required.",
              "\n2) `resistance` is deprecated. `data` should be used in its place.",
              "\n3) `p_mat` is deprecated. The P matrix is now provided via `data`.",
              "\n4) `tr_fun` is deprecated. See the samc() function documentation for the alternative.",
              "\n5) `directions` is deprecated. See the samc() function documentation for the alternative.",
              "\n6) `override` is deprecated. See the samc-class documentation for the alternative.",
              "\n\nPrevious code should continue to work except in an unlikely corner case, but version",
              "1.5.0 of the package will remove the deprecated functionality, so users",
              "should update previous work accordingly.")
  packageStartupMessage(msg)
}
