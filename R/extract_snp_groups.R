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
extract_snp_groups = function(pheno_gwas, endo_gwas, map, thresh, filter_hapmap, hapmap){

    ## Clean gwas column names and add unique id
    pheno_gwas_clean = format_gwas(pheno_gwas, map, filter_hapmap, hapmap)
    endo_gwas_clean = format_gwas(endo_gwas, map, filter_hapmap, hapmap)
    
    ## Filter gwas based on p-value
    pheno_gwas_thresh = pheno_gwas_clean[pheno_gwas_clean$P < thresh,]
    endo_gwas_thresh = endo_gwas_clean[endo_gwas_clean$P < thresh,]
    
    ## Create extract SNPs below p-value threshold for each trait
    df_pheno = data.frame(SNP = pheno_gwas_thresh$unique_id, trait = "pheno")
    df_endo = data.frame(SNP = endo_gwas_thresh$unique_id, trait = "endo")
    df_combined = rbind(df_pheno, df_endo)
    
    ## Create label as to whether SNP is pheno-only, endo-only, or both
    snps_assoc_tab <- table(df_combined$SNP)
    snps_assoc <- rep(NA_character_, length(snps_assoc_tab))
    names(snps_assoc) <- names(snps_assoc_tab)
    
    snps_assoc[snps_assoc_tab == 2] <- "both"
    snps_assoc[snps_assoc_tab == 1 & names(snps_assoc_tab) %in% df_pheno$SNP] <- "pheno_only"
    snps_assoc[snps_assoc_tab == 1 & names(snps_assoc_tab) %in% df_endo$SNP] <- "endo_only"
    
    # Drop any unassigned (shouldn't exist, but just in case)
    snps_assoc <- snps_assoc[!is.na(snps_assoc)]
    
    return(snps_assoc)

}