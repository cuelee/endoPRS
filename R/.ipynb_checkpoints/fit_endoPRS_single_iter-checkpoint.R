#' Function to apply the first step of the endoPRS method to generate a polygenic risk score model for a particular penalty.
#'
#' This function applies the first part of the endoPRS model. It runs a weighted lasso model that penalizes the SNPs differently based on
#' whether they are  associated with only the phenotype, only the endophenotype, or both. It fits only one set of penalties. This allows for
#' more efficient parallelization of the fit_endoPRS function. Must be used in combination with refit_endoPRS.
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
#' @param thresh An numeric vector corresponding to the p-value thresholds. This corresponds to the GWAS p-value threshold used to determine association with the phenotype and endophenotype.
#' @param grid An optional data frame of weights to use for the weighted penalty in the lasso model. The first column must be w2 and correspond to the weights applied to the group of SNPs associated with only the endophenotype and the second column must be w3 and correspond to the weights applied to the group of SNPs associated with both the phenotype and endophenotype. The set of SNPs associated with only the phenotype is given a weight of 1. If not provided, the original grid of weights from the endoPRS manuscript (0.1, 0.5, 1, 2, 10) will be used for both w2 and w3.
#' @param iteration A numerical value corresponding to which entry of the grid to fit. For example, 3 corresponds to fitting weights provided in the third row of the grid.
#' @param NCORES An optional value corresponding to the number of cores to use. Otherwise, the number of cores will be learned using nb_cores().
#' @param save_folder A path to a directory that files can be written to. The fitted model and its performance in the validation set will be saved to that directory.
#' @param save_models An optional boolean value. Indicates whether to save the intermediate models to disk. If set to true, the save_folder will be used.
#'
#' @return A list with two elements: the model itself and the performance in the validation set
#' \itemize{
#' \item{model: The final fitted endoPRS model as an object of class big_sp_list. This can be used to apply the final model to a test set using the predict() function from the bigsnpr package.  }
#' \item{val_perf: A data frame consisting of the performance of the model in the validation set. R^2 is used for linear models and AUC is used for logistic models. }
#' }
#'
#' @export
fit_endoPRS_single_iter = function(G, map, fam,
                       train_pheno, train_covar,
                       val_pheno, val_covar,
                       pheno_gwas, endo_gwas,
                       filter_hapmap = F, hapmap = NULL, type = NULL,
                       thresh,
                       grid = NULL, iteration, NCORES = NULL,
                       save_folder,
                       save_models = F){

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
  if(is.null(NCORES)){
    NCORES = nb_cores()
  }


  ## Create grid of weights
  if(is.null(grid)){
    w2 = c(1e-1, 0.5, 1, 2, 10)
    w3 = c(1e-1, 0.5, 1, 2,  10)
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
    hapmap$unique_id1 = paste(hapmap$CHR, hapmap$SNP, hapmap$A1, hapmap$A2, sep = "_")
    hapmap$unique_id2 = paste(hapmap$CHR, hapmap$SNP, hapmap$A2, hapmap$A1, sep = "_")
  }

  ## Add unique ID to map
  map$unique_id = paste(as.numeric(map$chromosome), map$marker.ID,
                        map$allele1, map$allele2, sep = "_")


  #######################################################################################
  ##################      Part 4: Fit Weighted Model over grid     ######################
  #######################################################################################


  print(paste("Fitting model corresponding to thresh:", thresh, "and penalties:" ))
  print(grid[iteration,])

  ## Obtain SNPs associated with pheno, with endo, and with both
  snps_assoc =  extract_snp_groups(pheno_gwas, endo_gwas, map, thresh, filter_hapmap, hapmap)

  ## Which SNPs to use
  ind.col = which(map$unique_id %in% names(snps_assoc))


  ## Create table of penalties to use
  penalty_table = create_penalty_table(snps_assoc, map, ind.col, w2 = grid$w2[iteration], w3 = grid$w3[iteration])

  ## Set seed for reproducibility
  set.seed(1)


  ## Fit endoPRS weighted model
  model = fit_weighted_model(G, y.train, train_index, ind.col, penalty_table, covar.train, NCORES, type)

  ## Save model if specified
  if(save_models){
    save_file = paste0(save_folder, "/trainingonly_model_thresh", thresh, "_w2", grid$w2[iteration],"_w3", grid$w3[iteration],".rds")
    save(model, file = save_file)
  }

  ## Apply to validation set
  pred_val <- predict(model, G, val_index, covar.row = covar.val)


  ## Performance in validation set
  if(type == "linear"){
    val_perf = cor(pred_val, y.val)^2
  }

  ## Performance in validation set
  if(type == "logistic"){
    val_perf = AUC(pred_val, y.val)
  }

  ## Save results
  res = data.frame(thresh = thresh, w2 = grid$w2[iteration],
                     w3 = grid$w3[iteration], val_res = val_perf)
  save_file = paste0(save_folder, "/trainingonly_model_thresh", thresh, "_w2", grid$w2[iteration],"_w3", grid$w3[iteration],"_validationresults.csv")
  write.csv(res, save_file, row.names = F)


  ## Return final model and final beta
  return(list(model = model , val_perf = res))

}
