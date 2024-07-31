
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ReducedExperiment

<!-- badges: start -->

[![R-CMD-check](https://github.com/jackgisby/ReducedExperiment/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jackgisby/ReducedExperiment/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/jackgisby/ReducedExperiment/graph/badge.svg?token=FHNH7AA6S3)](https://codecov.io/gh/jackgisby/ReducedExperiment)
[![GitHub
issues](https://img.shields.io/github/issues/jackgisby/ReducedExperiment)](https://github.com/jackgisby/ReducedExperiment/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/jackgisby/ReducedExperiment)](https://github.com/jackgisby/ReducedExperiment/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/ReducedExperiment.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/ReducedExperiment)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/ReducedExperiment.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/ReducedExperiment)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/ReducedExperiment.svg)](http://bioconductor.org/packages/stats/bioc/ReducedExperiment/)
[![Bioc
support](https://bioconductor.org/shields/posts/ReducedExperiment.svg)](https://support.bioconductor.org/tag/ReducedExperiment)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/ReducedExperiment.svg)](https://bioconductor.org/packages/release/bioc/html/ReducedExperiment.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/ReducedExperiment.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/ReducedExperiment/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/ReducedExperiment.svg)](https://bioconductor.org/packages/release/bioc/html/ReducedExperiment.html#since)
[![check-bioc](https://github.com/jackgisby/ReducedExperiment/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/jackgisby/ReducedExperiment/actions/workflows/check-bioc.yml)
<!-- badges: end -->

# ReducedExperiment <img src="inst/ReducedExperiment_hex.png" align="right" height="174" width="150" />

ReducedExperiment provides SummarizedExperiment-like containers for
storing and manipulating dimensionally-reduced assay data. The
ReducedExperiment classes allow users to simultaneously manipulate their
original dataset and their decomposed data, in addition to other
method-specific outputs like feature loadings. Implements utilities and
specialised classes for the application of stabilised independent
component analysis (sICA) and weighted gene correlation network analysis
(WGCNA).

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `ReducedExperiment`
from [Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("ReducedExperiment")
```

And the development version from
[GitHub](https://github.com/jackgisby/ReducedExperiment) with:

``` r
BiocManager::install("jackgisby/ReducedExperiment")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("ReducedExperiment")
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
#>     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     Position, rank, rbind, Reduce, rownames, sapply, setdiff, table,
#>     tapply, union, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> 
#> Attaching package: 'IRanges'
#> The following object is masked from 'package:grDevices':
#> 
#>     windows
#> Loading required package: GenomeInfoDb
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
#> 
#> Attaching package: 'ReducedExperiment'
#> The following objects are masked from 'package:Biobase':
#> 
#>     featureNames, featureNames<-, sampleNames, sampleNames<-
#> The following object is masked from 'package:stats':
#> 
#>     loadings
```

``` r
## basic example code
```

## Citation

Below is the citation output from using `citation('ReducedExperiment')`
in R. Please run this yourself to check for any updates on how to cite
**ReducedExperiment**.

``` r
print(citation("ReducedExperiment"), bibtex = TRUE)
#> To cite package 'ReducedExperiment' in publications use:
#> 
#>   jackgisby (2024). _ReducedExperiment_.
#>   doi:10.18129/B9.bioc.ReducedExperiment
#>   <https://doi.org/10.18129/B9.bioc.ReducedExperiment>,
#>   https://github.com/jackgisby/ReducedExperiment/ReducedExperiment - R
#>   package version 0.1.2,
#>   <http://www.bioconductor.org/packages/ReducedExperiment>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {ReducedExperiment},
#>     author = {{jackgisby}},
#>     year = {2024},
#>     url = {http://www.bioconductor.org/packages/ReducedExperiment},
#>     note = {https://github.com/jackgisby/ReducedExperiment/ReducedExperiment - R package version 0.1.2},
#>     doi = {10.18129/B9.bioc.ReducedExperiment},
#>   }
```

Please note that the `ReducedExperiment` was only made possible thanks
to many other R and bioinformatics software authors, which are cited
either in the vignette.

## Code of Conduct

Please note that the `ReducedExperiment` project is released with a
[Contributor Code of
Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.19/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation
  website](http://jackgisby.github.io/ReducedExperiment) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.19/biocthis)*.
