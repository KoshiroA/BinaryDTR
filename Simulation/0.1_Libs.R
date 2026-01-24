

if (requireNamespace("doParallel", quietly = TRUE)) {
  suppressWarnings(library(doParallel))
} else {
  stop("Required package 'doParallel' is not installed.")
}

if (requireNamespace("doRNG", quietly = TRUE)) {
  suppressWarnings(library(doRNG))
} else {
  stop("Required package 'doRNG' is not installed.")
}

if (requireNamespace("brm", quietly = TRUE)) {
  suppressWarnings(library(brm))
} else {
  stop("Required package 'brm' is not installed.")
}

if (requireNamespace("geeM", quietly = TRUE)) {
  suppressWarnings(library(geeM))
} else {
  stop("Required package 'geeM' is not installed.")
}

if (requireNamespace("tidyverse", quietly = TRUE)) {
  suppressWarnings(library(tidyverse))
} else {
  stop("Required package 'tidyverse' is not installed.")
}

if (requireNamespace("nleqslv", quietly = TRUE)) {
  suppressWarnings(library(nleqslv))
} else {
  install.packages("nleqslv")
}

if (requireNamespace("patchwork", quietly = TRUE)) {
  suppressWarnings(library(patchwork))
} else {
  install.packages("patchwork")
}

if (requireNamespace("mgcv", quietly = TRUE)) {
  suppressWarnings(library(mgcv))
} else {
  install.packages("mgcv")
}

if (requireNamespace("geepack", quietly = TRUE)) {
  suppressWarnings(library(geepack))
} else {
  stop("Required package 'geepack' is not installed.")
}

if (requireNamespace("boot", quietly = TRUE)) {
  suppressWarnings(library(boot))
} else {
  stop("Required package 'boot' is not installed.")
}

# library(ggplot2)
# library(gridExtra)
# library(openxlsx)
# library(patchwork)
