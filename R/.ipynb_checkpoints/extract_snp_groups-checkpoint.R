#' Function to extract snps associated with phenotype, with endotype, and with both
#'
#' This function creates a table that contains the name of the SNP and the corresponding
#' group it is assigned to based on GWAS p-values (pheno only, endo only, or both). It is
#' used to determine penalties for the weighted lasso model.
#'
#' @param pheno_gwas data frame containing results of phenotype GWAS. It should contain a column with chromosome, SNP ID, allele1, allele2, and P-value.
#' @param endo_gwas data frame containing results of endophenotype GWAS. It should contain a column with chromosome, SNP ID, allele1, allele2, and P-value.
#' @param map data frame containing information about SNPs in G (genotype matrix)
#' @param thresh P-value threshold used to determine SNP association with phenotype and endophenotype in GWAS.
#' @param filter_hapmap logical value as to whether to filter to SNPs in hapmap.
#' @param hapmap data frame with hapmap variants. Only provide if filter_hapmap argument set to TRUE. It should contain columns CHR, SNP, A1, A2.
#'
#' @return A named table. The name corresponds to the SNP and value corresponds to the assigned group (pheno only, endo only, or both). It is used as input in the create_penalty_table function

## Function to extract snps associated with pheno, with endo, and with both
extract_snp_groups <- function(pheno_gwas, endo_gwas, map, thresh, filter_hapmap, hapmap) {
  # Clean and deduplicate GWAS inputs
  pheno_gwas_clean <- format_gwas(pheno_gwas, map, filter_hapmap, hapmap)
  endo_gwas_clean  <- format_gwas(endo_gwas,  map, filter_hapmap, hapmap)

  pheno_gwas_clean <- pheno_gwas_clean[!duplicated(pheno_gwas_clean$unique_id), ]
  endo_gwas_clean  <- endo_gwas_clean[!duplicated(endo_gwas_clean$unique_id),  ]

  # Extract significant SNPs
  pheno_sig <- pheno_gwas_clean$unique_id[pheno_gwas_clean$P < thresh]
  endo_sig  <- endo_gwas_clean$unique_id[endo_gwas_clean$P < thresh]

  # Create named character vector efficiently
  snps_assoc <- character()
  snps_assoc[pheno_sig] <- "pheno_only"
  snps_assoc[endo_sig]  <- ifelse(endo_sig %in% names(snps_assoc), "both", "endo_only")

  # Sanity check
  allowed <- c("pheno_only", "endo_only", "both")
  if (any(!snps_assoc %in% allowed)) stop("Unexpected values in snps_assoc.")

  return(snps_assoc)
}

