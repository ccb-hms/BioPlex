# BioPlex
Analysis of PPI data from Gygi lab

## Installation

Make sure to have the latest release version of 
[R](https://cran.r-project.org/) and 
[Bioconductor](https://bioconductor.org/install/) installed.

Then proceed from within R via:

```
BiocManager::install("ccb-hms/BioPlex", 
                     build_vignettes = TRUE,
                     auth_token = <your_auth_token>)
```

NOTE: you will need the `remotes` package to install from github. 
You will also need to provide the `auth_token` argument to 
`remotes::install_github` as the repo is private.        

Once you have the package installed, you can inspect the vignettes via:

```
browseVignettes("BioPlex")
```

## CORUM

See [here](https://github.com/ccb-hms/BioPlex/blob/0ca36e34957a4e7b0d34ee66915e5f4e5989cee4/vignettes/BioPlex.Rmd#L16) for how to access the CORUM protein complex data from within the package.

