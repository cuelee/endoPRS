library(purrr)
library(bigsnpr)
library(data.table)
library(tidyverse)
library(endoPRS)


#########################################################################################
#################################   Load parameters  ####################################
#########################################################################################

folder = commandArgs(trailingOnly=TRUE)[1] # folder to save results
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
data(pheno_gwas_refit) # GWAS of phenotype ran on combined training and validation set
data(endo_gwas_refit) # GWAS of endophenotype ran on combined training and validation set

## Load phenotype information of training, validation, and testing individuals
data(train_pheno) # Training set phenotypes
data(val_pheno) # Validation set phenotypes

## Load covariate information of training, validation, and testing individuals
data(train_covar) # Training set covariates
data(val_covar) # Validation set covariates


#########################################################################################
###################################   Refit Model    ######################################
#########################################################################################

refit_endoPRS(G, map, fam,
                         train_pheno, train_covar,
                         val_pheno, val_covar,
                         pheno_gwas = pheno_gwas_refit, endo_gwas = endo_gwas_refit,
                         filter_hapmap = F, hapmap = NULL, type = NULL,
                         threshes = c(0.01, 1e-4, 1e-6), NCORES = 5,
                         save_folder = folder,
                         save_model = T)
