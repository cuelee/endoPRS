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
  
  ## Clean GWAS column names and add unique ID
  pheno_gwas_clean <- format_gwas(pheno_gwas, map, filter_hapmap, hapmap)
  endo_gwas_clean  <- format_gwas(endo_gwas,  map, filter_hapmap, hapmap)
  
  ## Remove duplicated SNPs based on (CHR, BP)
  pheno_gwas_clean <- pheno_gwas_clean[!duplicated(pheno_gwas_clean[, c("CHR", "BP")]), ]
  endo_gwas_clean  <- endo_gwas_clean[!duplicated(endo_gwas_clean[, c("CHR", "BP")]), ]
  
  ## Filter GWAS by p-value threshold
  pheno_gwas_thresh <- pheno_gwas_clean[pheno_gwas_clean$P < thresh, ]
  endo_gwas_thresh  <- endo_gwas_clean[endo_gwas_clean$P < thresh, ]
  
  ## Label SNPs by trait
  df_pheno    <- data.frame(SNP = pheno_gwas_thresh$unique_id, trait = "pheno")
  df_endo     <- data.frame(SNP = endo_gwas_thresh$unique_id, trait = "endo")
  df_combined <- rbind(df_pheno, df_endo)
  
  ## Assign SNP groups
  snps_assoc_tab <- table(df_combined$SNP)
  snps_assoc     <- rep(NA_character_, length(snps_assoc_tab))
  names(snps_assoc) <- names(snps_assoc_tab)
  
  snps_assoc[snps_assoc_tab == 2] <- "both"
  snps_assoc[snps_assoc_tab == 1 & names(snps_assoc_tab) %in% df_pheno$SNP] <- "pheno_only"
  snps_assoc[snps_assoc_tab == 1 & names(snps_assoc_tab) %in% df_endo$SNP]  <- "endo_only"
  
  ## Drop any unassigned (should not exist, but safe)
  snps_assoc <- snps_assoc[!is.na(snps_assoc)]
  
  ## Sanity check
  if (!all(snps_assoc %in% c("pheno_only", "endo_only", "both"))) {
    stop("snps_assoc contains invalid group labels.")
  }
  
  return(snps_assoc)
}
