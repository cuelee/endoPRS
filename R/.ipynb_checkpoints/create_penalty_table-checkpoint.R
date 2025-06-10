#' Function to create table of penalties
#'
#' This function creates a data frame that contains the name of the SNP,
#' the group (pheno only, endo only, or both), and the associated penalty
#' to be used in the weighted lasso model.
#'
#' @param snps_assoc output of the extract_snp_groups function: a named table where the names are SNPs and values consist of whether SNP is pheno only, endo only, or associated with both.
#' @param map data frame containing information about SNPs in G (genotype matrix)
#' @param ind.col index of SNPs in map that are in snp_assoc table
#' @param w2 weight to assign SNPs in endo only group
#' @param w3 weight to assign SNPs associated with both pheno and endo
#'
#' @return a data frame that contains the name of the SNP, the group (pheno only, endo only, or both), and the associated penalty to be used in the weighted lasso model.

## Function to create table of penalties
create_penalty_table = function(snps_assoc, map, ind.col, w2, w3){

  ## Create penalty table - with SNP IDs and with group it is
  penalty_table = data.frame(RSID = names(snps_assoc)[match(map$unique_id[ind.col],names(snps_assoc))],
                             group = as.character(snps_assoc)[match(map$unique_id[ind.col],names(snps_assoc))])

  ## Create penalty factor - w2 for endo only, and w3 for both
  penalty_table$penalty = NA
  penalty_table$penalty[penalty_table$group == "pheno_only"] = 1
  penalty_table$penalty[penalty_table$group == "endo_only"] = w2
  penalty_table$penalty[penalty_table$group == "both"] = w3

  ## Return penalty table
  return(penalty_table)
}

