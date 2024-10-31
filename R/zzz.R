# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

# .onLoad <- function(libname, pkgname) {
#   assign("samc_env", new.env(parent = emptyenv()), parent.env())
# }

.onAttach <- function(libname, pkgname) {
  msg = paste("Old code may be affected by breaking changes. Check the package website and release notes for details.")
  packageStartupMessage(msg)
}
