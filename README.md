# endoPRS
R package to fit the endoPRS method. endoPRS is a weighted penalized regression model that incorporates information from endophenotypes to improve polygenic risk score prediction.


## Installation
To install endoPRS, you can use the following code:
```
library(devtools)
devtools::install_github("https://github.com/ekhar17/endoPRS")
library(endoPRS)
```

## Usage
The main function to fit endoPRS is the `fit_endoPRS` function. An example on how to run it:
```
endoPRS <- fit_endoPRS(G, map, fam, 
                       train_pheno, train_covar,
                       val_pheno, val_covar,
                       pheno_gwas, endo_gwas,
                       filter_hapmap = F,  hapmap = NULL, 
                       pheno_gwas_refit, endo_gwas_refit, 
                       save_folder = NULL)
```
Further information can be found in the vignette. 


## Citations
Please cite:
Kharitonova, E.V., et. al. EndoPRS: Incorporating Endophenotype Information to Improve Polygenic Risk Scores for Clinical Endpoints. *medRxiv*. (2024).
