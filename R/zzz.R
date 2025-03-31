.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    packageStartupMessage(
      "tvcure - Copyright (C) 2025 - Philippe LAMBERT\n",
      "This program comes with ABSOLUTELY NO WARRANTY; for details type `show_w()'.\n",
      "This is free software, and you are welcome to redistribute it\n",
      "under certain conditions; type `show_c()' for details."
    )
  }
}

#' Show Copyright Information
#'
#' Displays the copyright information for the package in interactive mode.
#'
#' @export
show_c <- function() {
  cat("This program is free software; you can redistribute it under certain conditions.\n")
  cat("For details, see the General Public License (GPL-3).\n")
}

#' Show Warranty Disclaimer
#'
#' Displays a notice that the program comes with absolutely no warranty.
#'
#' @export
show_w <- function() {
  cat("This program comes with ABSOLUTELY NO WARRANTY.\n")
  cat("For details, see the General Public License (GPL-3).\n")
}
