# BioPlex
R-side access to PPI data from Gygi lab

## Citation

If you use the BioPlex R package in published research, please cite:

Ludwig Geistlinger, Roger Vargas Jr, Tyrone Lee, Joshua Pan, Edward Huttlin, Robert Gentleman (2023) BioPlexR and BioPlexPy: integrated data products for the analysis of human protein interactions. *Bioinformatics*. doi: [10.1093/bioinformatics/btad091](https://doi.org/10.1093/bioinformatics/btad091).

## Installation

### Installation from Bioconductor

Users interested in using the stable release version of the `BioPlex` 
package: please follow the installation instructions
[here](https://bioconductor.org/packages/BioPlex). 
This is the recommended way of installing the package.

### Installation from GitHub

It is also possible to install the package directly from GitHub. This is 
for users/developers interested in using the latest development version of the
`BioPlex` package. 

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
