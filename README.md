# BioPlex
R-side access to PPI data from Gygi lab

## Installation

Make sure to have the latest release version of 
[R](https://cran.r-project.org/) and 
[Bioconductor](https://bioconductor.org/install/) installed.

Then proceed from within R via:

```
BiocManager::install("ccb-hms/BioPlex", 
                     auth_token = <your_auth_token>)
```

NOTE: you will need the `remotes` package to install from github. 
You will also need to provide the `auth_token` argument to 
`remotes::install_github` as the repo is private.        

To build the package vignettes upon installation use:

```
BiocManager::install("ccb-hms/BioPlex",
                     build_vignettes = TRUE,
                     dependencies = TRUE,
                     auth_token = <your_auth_token>)
```

Once you have the package installed, you can inspect the vignettes from within
R via:

```
browseVignettes("BioPlex")
```


## Bioplex PPIs

See [here](https://github.com/ccb-hms/BioPlex/blob/1eb7c039045a6d710bde7c8017cccedb98d5e9b5/vignettes/BioPlex.Rmd#L29) for how to access the BioPlex protein-protein interaction data from within the package.

## CORUM complexes

See [here](https://github.com/ccb-hms/BioPlex/blob/0ca36e34957a4e7b0d34ee66915e5f4e5989cee4/vignettes/BioPlex.Rmd#L16) for how to access the CORUM protein complex data from within the package.

## 293T transcriptome data

### GSE122425 

#### Alignment with STAR

O2 project directory for code and data:

`/n/shared_db/ccb/bioplex`

The fastq files for SRP168405 were aligned with STAR and the following human genome:

`/n/groups/shared_databases/star_reference/hg19.genes.gtf
/n/groups/shared_databases/genomes/hg19.fa`

Alignment bam files are here: 

`/n/shared_db/ccb/bioplex/data/2_alignment/bam`

Samples GSM3466389, GSM3466390, GSM3466391 are wild-type, while samples GSM3466392, GSM3466393 and GSM3466394 are the NSUN2 knockouts.

STAR was run with default thresholds. Input parameters can be viewed here: 

`/n/shared_db/ccb/bioplex/code/2_alignment/2_align.sh`

The STAR genome index is saved in our group folder: 

`/n/shared_db/ccb/ref/star_indices/hg19`

Each sample's fastq file was ~30 GB, and each alignment took ~15 minutes of wall time and 30 GB of RAM, running with 16 threads on O2.

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
