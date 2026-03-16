# Contributing

Thank you for your interest in contributing to `integrative-drugrep-tb`!

## Setting up the environment

This project uses [renv](https://rstudio.github.io/renv/) to manage R dependencies.

```r
# Install renv if needed
install.packages("renv")

# Restore the project library
renv::restore()
```

R 4.4.2 and Bioconductor 3.20 are required. The Python notebook (`scripts/00_cleanup_expression_data.ipynb`) requires Python 3 with standard scientific libraries.

## Code style

- R code should be formatted with [styler](https://styler.r-lib.org/) and checked with [lintr](https://lintr.r-lib.org/) before submitting a PR.
- Function documentation follows [roxygen2](https://roxygen2.r-lib.org/) style (`#' @description`, `#' @param`, `#' @returns`).

## Submitting changes

1. Fork the repository and create a branch from `main`.
2. Make your changes and ensure the relevant scripts run end-to-end.
3. Open a pull request using the provided template and tag a reviewer.

## Questions

Open a [GitHub Issue](https://github.com/JRaviLab/integrative-drugrep-tb/issues) for bugs or feature requests, or contact the [JRaviLab](https://jravilab.github.io/) or [KrishnanLab](https://www.thekrishnanlab.org).
