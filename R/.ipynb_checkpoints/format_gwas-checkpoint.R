#' Function to clean and format GWAS results
#'
#' This function cleans a dataframe of GWAS results. It formats the column names to be readable by later functions. It adds
#' a unique id (combination of chromosome, SNP ID, allele1, and allele2). It keeps only the SNPs present in the
#' genotype data (those in the map data frame). And it can also filter to keep only hapmap variants.
#'
#' @param df data frame containing results of GWAS. It should contain a column with chromosome, SNP, allele1, allele2, and P-value.
#' @param map data frame containing information about SNPs in G (genotype matrix)
#' @param filter_hapmap logical value as to whether to filter to SNPs in hapmap.
#' @param hapmap data frame with hapmap variants. Only provide if filter_hapmap argument set to TRUE. It should contain columns CHR, SNP, A1, A2.
#'
#' @return A data frame consisting of GWAS results that include columns CHR, SNP, A0, A1, P and a unique id. It only contains SNPs that are present in the genotype data.

## Function to change column names of GWAS and add unique id 1 and 2 column
format_gwas = function(df, map, filter_hapmap, hapmap){
    ## Properly format chromosome, snp, a0, a1 and p-value column for ease of use
    colnames(df)[toupper(colnames(df) ) %in% c("CHROM", "CHR")] = "CHR"
    colnames(df)[toupper(colnames(df) ) %in% c("ID", "RSID", "SNP")] = "SNP"
    colnames(df)[toupper(colnames(df) ) %in% c("ALLELE0", "A0", "ALLELE2", "A2")] = "A0"
    colnames(df)[toupper(colnames(df) ) %in% c("ALLELE1", "A1")] = "A1"
    colnames(df)[toupper(colnames(df) ) %in% c("P", "PVAL")] = "P"
    
    ## Create unique ID column
    unique_id1 = paste(df$CHR, df$BP, df$A0, df$A1, sep = "_")
    unique_id2 = paste(df$CHR, df$BP, df$A1, df$A0, sep = "_")
    
    ## Add the version that matches map
    df$unique_id = unique_id1
    df$unique_id[unique_id2 %in% map$unique_id] = unique_id2[unique_id2 %in% map$unique_id]
    
    ## Keep only those snps that overlap with map
    df = df[df$unique_id %in% map$unique_id,]
    
    ## Create unique ID column
    unique_id1_hapmap = paste(df$CHR, df$BP, sep = "_")
    unique_id2_hapmap = paste(df$CHR, df$BP, sep = "_")
    df$unique_id_hapmap = unique_id1_hapmap
    
    ## If filter hapmap options selected, keep only those that overlap with hapmap
    if(filter_hapmap){
    df = df[c(df$unique_id_hapmap %in% hapmap$unique_id),]
    }

    ## Return data frame with proper column names
    return(df)
}
