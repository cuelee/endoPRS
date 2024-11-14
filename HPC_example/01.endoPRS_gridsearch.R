library(purrr)
library(bigsnpr)
library(data.table)
library(tidyverse)
#library(endoPRS)


#########################################################################################
#################################   Load parameters  ####################################
#########################################################################################

thresh = commandArgs(trailingOnly=TRUE)[1] # p-value threshold to use
iteration = commandArgs(trailingOnly=TRUE)[2] # iteration of grid to use
folder = commandArgs(trailingOnly=TRUE)[3] # folder to save results
thresh = as.numeric(thresh)
iteration = as.numeric(iteration)
print(paste("Thresh is",thresh))
print(paste("Iteration is",iteration))
print(paste("Folder to save results is",folder))


#########################################################################################
###################################   Load Data    ######################################
#########################################################################################

## Load bigsnpr data and load the example genotype data
bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
rds <- snp_readBed(bedfile, backingfile = tempfile())
test <- snp_attach(rds)

## Specify the genotype, snp information, and participant information
G = test$genotypes # Genotypes
map = test$map # SNP information
fam  = test$fam # Participant information

## Load GWAS information for phenotype and endophenotype
data(pheno_gwas) # GWAS of phenotype ran on training set
data(endo_gwas) # GWAS of endophenotype ran on training set

## Load phenotype information of training, validation, and testing individuals
data(train_pheno) # Training set phenotypes
data(val_pheno) # Validation set phenotypes

## Load covariate information of training, validation, and testing individuals
data(train_covar) # Training set covariates
data(val_covar) # Validation set covariates



#########################################################################################
###################################   Fit Model    ######################################
#########################################################################################

fit_endoPRS_single_iter(G, map, fam,
                       train_pheno, train_covar,
                       val_pheno, val_covar,
                       pheno_gwas, endo_gwas,
                       filter_hapmap = F, hapmap = NULL, type = NULL,
                       thresh,
                       grid = NULL, iteration,
                       NCORES = 5,
                       save_folder = folder,
                       save_models = F)
