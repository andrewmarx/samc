# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

# .onLoad <- function(libname, pkgname) {
#   assign("samc_env", new.env(parent = emptyenv()), parent.env())
# }

.onAttach <- function(libname, pkgname) {
  msg <-paste("Version 2 and Version 3 had breaking changes. Check the package website for details.")
  packageStartupMessage(msg)
}
