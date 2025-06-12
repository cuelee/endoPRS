#' Function to apply endoPRS method to generate a polygenic risk score model.
#'
#' This function applies the endoPRS model. It runsa weighted lasso model that penalizes the SNPs differently based on
#' whether they are  associated with only the phenotype, only the endophenotype, or both. The optimal set of weights is
#' determined by validation set performance. The model is then refit with the chosen weights on the combined training and
#' validation set to obtain the final PRS model.
#'
#'
#' @param G An object of class FBM from the bigsnpr package that contains genotypes of the individuals used for training and those used for validation.
#' @param map A data frame containing information about SNPs in G (genotype matrix). This can be the map object form the bigSNP class. It should contain columns labeled chromosome, marker ID, and allele1 and allele2.
#' @param fam A data frame containing information about the individuals in G (genotype matrix). Column one must correspond to FID and column 2 must correspond to IID.
#' @param train_pheno A data frame with three columns corresponding to the phenotypes of individuals used for training. Column one is FID, column two is IID, and column three is y the phenotype. If the phenotype is binary, y must consist of 0's and 1's.
#' @param train_covar A data frame corresponding with the covariates of the individuals used for training (ie genetic PC's, sex, age, etc). The first column must be FID, and the second column must be IID. The order of individuals must be the same as train_pheno. Covariates that are factors (i.e. assessment center) must be encoded as such. The covariates are not penalized in the endoPRS model.
#' @param val_pheno A data frame with three columns corresponding to the phenotypes of individuals used for validation. Column one is FID, column two is IID, and column three is y the phenotype. If the phenotype is binary, y must consist of 0's and 1's.
#' @param val_covar A data frame corresponding with the covariates of the individuals used for validation (ie genetic PC's, sex, age, etc). The first column must be FID, and the second column must be IID. The order of individuals must be the same as val_pheno. Covariates that are factors (i.e. assessment center) must be encoded as such. The covariates are not penalized in the endoPRS model.
#' @param pheno_gwas A data frame containing the results of the GWAS run on the phenotype. It must not include any individuals from the validation set. It must contain columns corresponding to chromosome, SNP, allele1, allele2, and P-value.
#' @param endo_gwas A data frame containing the results of the GWAS run on the endophenotype. It must not include any individuals from the validation set. It must contain columns corresponding to chromosome, SNP, allele1, allele2, and P-value.
#' @param filter_hapmap An optional logical character. It corresponds to whether to only run the model on variants in the hapmap3 set.
#' @param hapmap An optional data frame with hapmap variants. This must be provided if filter_hapmap is set to TRUE. The hapmap3 data frame should contain columns CHR, SNP, A1, A2.
#' @param type An optional character vector of "linear" or "logistic." It specifies which type of model to fit to the data. If not provided, it will be learned from the phenotype. If the phenotype consists of 0's and 1's penalized logistic regression will be used, otherwie penalized linear regression will be used.
#' @param threshes An optional character of p-value thresholds. This corresponds to the GWAS p-value thresholds used to determine association with the phenotype and endophenotype. The optimal p-value threshold is learned from tuning on the validation set. If not provided, the original values from the endoPRS manuscript will be used (0.01, 1e-4, and 1e-6).
#' @param grid An optional data frame of weights to use for the weighted penalty in the lasso model. The first column must be w2 and correspond to the weights applied to the group of SNPs associated with only the endophenotype and the second column must be w3 and correspond to the weights applied to the group of SNPs associated with both the phenotype and endophenotype. The set of SNPs associated with only the phenotype is given a weight of 1. If not provided, the original grid of weights from the endoPRS manuscript (0.1, 0.5, 1, 2, 10) will be used for both w2 and w3.
#' @param NCORES An optional value corresponding to the number of cores to use. Otherwise, the number of cores will be learned using nb_cores().
#' @param pheno_gwas_refit An optional data frame containing the results of the GWAS run on the phenotype using both the training and validation set. It must contain columns corresponding to chromosome, SNP, allele1, allele2, and P-value.
#' @param endo_gwas_refit An optional data frame containing the results of the GWAS run on the endophenotype using both the training and validation set. It must contain columns corresponding to chromosome, SNP, allele1, allele2, and P-value.
#' @param save_folder An optional argument. If included, it should specify a path to a directory that files can be written to. If specified, all the models that are fit are saved to that directory. This is particularly useful for larger data sets, as model fitting can take a while, so saving intermediate results can prevent the need for rerunning the same models again.
#'
#' @return A list with two elements: the betas of the final model and the model itself
#' \itemize{
#' \item{beta: A data frame consisting of the SNPs included in the final model and their corresponding coefficients. This can be used to apply the PRS to an external data set using software such as PLINK.}
#' \item{model: The final fitted endoPRS model as an object of class big_sp_list. This can be used to apply the final model to a test set using the predict() function from the bigsnpr package.  }
#' }
#'
#' @export
fit_endoPRS = function(G, map, fam,
                       train_pheno, train_covar,
                       val_pheno, val_covar,
                       pheno_gwas, endo_gwas,
                       filter_hapmap = F, hapmap = NULL, type = NULL,
                       threshes = c(1e-2, 1e-4, 1e-6),
                       grid = NULL, NCORES = NULL,
                       pheno_gwas_refit, endo_gwas_refit,
                       save_folder = NULL){

    bigstatsr::assert_cores(NCORES)  # Enforce safe core usage
    
    ## Check that there are no missing values
    if(any(is.na(train_pheno)) |  any(is.na(val_pheno)) |
     any(is.na(train_covar)) | any(is.na(val_covar)) ){
    stop("Currently endoPRS cannot handle any missing values in the training or validation phenotype or covariates. Please reformat.")
    }
    
    ## Check that ID's of train pheno & covar match, as well as that for validation set
    if(!all.equal(train_pheno[,1:2], train_covar[,1:2])){
    stop("First two columns of train pheno and train covar do not match. Please reformat.")
    }
    if(!all.equal(val_pheno[,1:2], val_covar[,1:2])){
    stop("First two columns of val pheno and val covar do not match. Please reformat.")
    }
    
    ## Check G, fam, and map input
    if(nrow(G) != nrow(fam) | ncol(G) != nrow(map)){
    stop("Number of rows of G not match the number of rows of fam. Please check your input.")
    }
    if(ncol(G) != nrow(map)){
    stop("Number of columns of G not match the number of rows of map. Please check your input.")
    }
    
    ## Parallelize
    if (is.null(NCORES)) {
      NCORES <- bigstatsr::nb_cores()
    }
    
    ## Create grid of weights
    if(is.null(grid)){
    w2 = c(1e-1, 0.5, 1, 2, 10)
    w3 = c(1e-1, 0.5, 1, 2, 10)
    grid = expand.grid(w2 = w2, w3 = w3)
    }
    
    #######################################################################################
    ###########     Part 2: Process Phenotype and Covariates    ###########################
    #######################################################################################
    
    
    ## Obtain ID's of all individuals, training set, and validation set
    geno_id = paste(fam[,1], fam[,2], sep = "_")
    train_id = paste(train_pheno[,1], train_pheno[,2], sep = "_")
    val_id = paste(val_pheno[,1], val_pheno[,2], sep = "_")
    
    
    ## Check that there is overlap between ID's
    if( sum(train_id %in% geno_id) == 0 | sum(val_id %in% geno_id) == 0  ){
    stop("There is no overlap between training/validation ID's and genotype ID's. Please reformat.")
    }
    
    ## Obtain indices in geno id that correspond to train and validation
    train_index = match(train_id, geno_id)
    val_index = match(val_id, geno_id)
    
    ## Obtain y-train and y-val
    y.train = c(train_pheno[,3])
    y.val = c(val_pheno[,3])
    
    ## Infer type or regression if not provided
    if(is.null(type)){
        if(sum(!(y.train %in% c(0,1))) > 0){
          type = "linear"
        } else {
          type = "logistic"
        }
    }
    
    ## Obtain Covariates for training and validation
    covar.train =covar_from_df(train_covar[,-c(1,2)])
    covar.val =covar_from_df(val_covar[,-c(1,2)])    
    
    
    #######################################################################################
    ######################   Part 3: Process Variant Info   ###############################
    #######################################################################################
    
    
    ## Filter to hapmap variants
    if(filter_hapmap){
    if(is.null(hapmap)){
      stop("Must provide hapmap data frame if selected filter_hapmap option.")
    }
    hapmap$unique_id = paste(hapmap$CHR, hapmap$BP, sep = "_")
    }
    
    ## Add unique ID to map
    map$unique_id = paste(as.numeric(map$chromosome), map$physical.pos,
                        map$allele1, map$allele2, sep = "_")

    
    #######################################################################################
    ##################      Part 4: Fit Weighted Model over Grid     ######################
    #######################################################################################
    
    val_results <- data.frame()
    
    cat("Fitting grid of", nrow(grid), "models over threshes:", paste(threshes, collapse = ", "), "\n")
    
    for (thresh in threshes) {
    
      # Extract SNP associations
      snps_assoc <- extract_snp_groups(pheno_gwas, endo_gwas, map, thresh, filter_hapmap, hapmap)
    
      # Filter and match SNPs
      valid_ids  <- intersect(map$unique_id, names(snps_assoc))
      snps_assoc <- snps_assoc[valid_ids]
      ind.col    <- which(map$unique_id %in% valid_ids)
    
      cat("Using", length(ind.col), "SNPs for threshold", thresh, "\n")
    
      pb <- txtProgressBar(min = 0, max = nrow(grid), style = 3)
      for (iter in seq_len(nrow(grid))) {
    
        # Generate penalty table
        penalty_table <- create_penalty_table(snps_assoc, map, ind.col, w2 = grid$w2[iter], w3 = grid$w3[iter])
    
        set.seed(1)
        model <- fit_weighted_model(G, y.train, train_index, ind.col, penalty_table, covar.train, NCORES, type)
    
        # Save model if path specified
        if (!is.null(save_folder)) {
          save_file <- file.path(save_folder, sprintf("trainingonly_model_thresh%s_w2%s_w3%s.rds", thresh, grid$w2[iter], grid$w3[iter]))
          save(model, file = save_file)
        }
    
        # Predict and evaluate
        pred_val <- predict(model, G, val_index, covar.row = covar.val)
    
        val_perf <- if (type == "linear") {
          cor(pred_val, y.val)^2
        } else if (type == "logistic") {
          AUC(pred_val, y.val)
        } else {
          stop("Invalid model type: must be 'linear' or 'logistic'")
        }
    
        val_results <- rbind(val_results, data.frame(thresh = thresh, w2 = grid$w2[iter],
                                                     w3 = grid$w3[iter], val_res = val_perf))
    
        setTxtProgressBar(pb, iter)
      }
    
      close(pb)
      cat("Finished threshold:", thresh, "\n")
    }
    
    if (!is.null(save_folder)) {
      write.csv(val_results, file = file.path(save_folder, "trainingonly_models_results.csv"), row.names = FALSE)
    }


    #######################################################################################
    ###########     Part 5: Refit using combined training + validation     ################
    #######################################################################################
    
    ## Determine which model is best performing
    best_performing_model  = val_results[which.max(val_results$val_res),]
    
    print(paste("Best Peforming Model Corresponds to Thresh:", best_performing_model$thresh,
              "w2:", best_performing_model$w2, "w3:", best_performing_model$w3 ))
    print("Start refitting model with combined training and validation set.")
    
    ## If alternate gwasval_results using combined training+validation set provided, use those SNPs for refitting
    if(!is.null(pheno_gwas_refit) & !is.null(endo_gwas_refit) ){
        pheno_gwas = pheno_gwas_refit
        endo_gwas = endo_gwas_refit
        print("Using refit summary statistics")
    }
    
    ## Set thresh and weights for refitting
    thresh = best_performing_model$thresh
    w2 = best_performing_model$w2
    w3 = best_performing_model$w3
    
    ## Obtain SNPs associated with pheno, with endo, and with both
    snps_assoc =  extract_snp_groups(pheno_gwas, endo_gwas, map, thresh, filter_hapmap, hapmap)
    
    ## Which SNPs to use
    ind.col = which(map$unique_id %in% names(snps_assoc))
    
    ## Create table of penalties to use
    penalty_table = create_penalty_table(snps_assoc, map, ind.col, w2 = w2, w3 = w3)
    
    ## Refit endoPRS weighted model
    model = fit_weighted_model(G, c(y.train,y.val), c(train_index,val_index),
                             ind.col, penalty_table, rbind(covar.train, covar.val), NCORES, type)
    ## If save folder specified, save model
    if(!is.null(save_folder)){
    save_file = paste0(save_folder, "/combined_trainval_thresh", thresh, "_w2", w2,"_w3", w3,".rds")
    save(model, file = save_file)
    }
    
    ## Extract beta coefficients from fitted model
    beta_info = obtain_beta(model, ind.col, map)
    
    ## Return final model and final beta
    return(list(beta = beta_info, model = model))

}