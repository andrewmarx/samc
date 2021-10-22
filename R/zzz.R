
# .onLoad <- function(libname, pkgname) {
#   assign("samc_env", new.env(parent = emptyenv()), parent.env())
# }

.onAttach <- function(libname, pkgname) {
  msg <-paste("The following parameters have officially been removed from the samc() function in v1.5.0:",
              "\n1) `latlon`: Now handled automatically.",
              "\n2) `resistance`: `data` should be used in its place.",
              "\n3) `p_mat`: The P matrix is now provided via `data`.",
              "\n4) `tr_fun`: See the samc() function documentation for the alternative.",
              "\n5) `directions`: See the samc() function documentation for the alternative.",
              "\n6) `override`: See the samc-class documentation for the alternative.")
  packageStartupMessage(msg)
}
