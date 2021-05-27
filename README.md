# BioPlex
R-side access to PPI data from Gygi lab

## Installation

Make sure to have the latest release version of 
[R](https://cran.r-project.org/) and 
[Bioconductor](https://bioconductor.org/install/) installed.

Then proceed from within R via:

```
BiocManager::install("ccb-hms/BioPlex") 
```

NOTE: you will need the `remotes` package to install from github. 

To build the package vignettes upon installation use:

```
BiocManager::install("ccb-hms/BioPlex",
                     build_vignettes = TRUE,
                     dependencies = TRUE)
```

Once you have the package installed, you can inspect the vignettes from within
R via:

```
browseVignettes("BioPlex")
```

## Cleaning your cache 

Note that calling functions like `getCorum` or `getBioPlex` with argument
`cache = FALSE` will automatically overwrite the corresponding object in your 
cache. It is thus typically not required for a user to interact with the cache.

For more extended control of the cache, use from within R:

```
cache.dir <- tools::R_user_dir("BioPlex", which = "cache") 
bfc <- BiocFileCache::BiocFileCache(cache.dir)
```

and then proceed as described in the
[BiocFileCache vignette, Section 1.10](https://www.bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html#cleaning-or-removing-cache)

either via `cleanbfc()` to clean or `removebfc()` to remove your cache.

To do a hard reset (use with caution!):

```
BiocFileCache::removebfc(bfc)
```
