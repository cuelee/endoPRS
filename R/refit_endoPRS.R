#' Function to apply the second (refitting) step of the endoPRS method to generate a polygenic risk score model.
#'
#' This function applies the refitting step of the endoPRS model. It loads the results of fit_endoPRS_single_iter to determine
#' the optimal chosen weights and GWAS p-value threshold in the validation set. It then refits the endoPRS model using these tuning
#' parameters on the combined training and validation set to obtain the final PRS model. Must be used in combination with fit_endoPRS_single_iter.
#'
#'
#' @param G An object of class FBM from the bigsnpr package that contains genotypes of the individuals used for training and those used for validation.
#' @param map A data frame containing information about SNPs in G (genotype matrix). This can be the map object form the bigSNP class. It should contain columns labeled chromosome, marker ID, and allele1 and allele2.
#' @param fam A data frame containing information about the individuals in G (genotype matrix). Column one must correspond to FID and column 2 must correspond to IID.
#' @param train_pheno A data frame with three columns corresponding to the phenotypes of individuals used for training. Column one is FID, column two is IID, and column three is y the phenotype. If the phenotype is binary, y must consist of 0's and 1's.
#' @param train_covar A data frame corresponding with the covariates of the individuals used for training (ie genetic PC's, sex, age, etc). The first column must be FID, and the second column must be IID. The order of individuals must be the same as train_pheno. Covariates that are factors (i.e. assessment center) must be encoded as such. The covariates are not penalized in the endoPRS model.
#' @param val_pheno A data frame with three columns corresponding to the phenotypes of individuals used for validation. Column one is FID, column two is IID, and column three is y the phenotype. If the phenotype is binary, y must consist of 0's and 1's.
#' @param val_covar A data frame corresponding with the covariates of the individuals used for validation (ie genetic PC's, sex, age, etc). The first column must be FID, and the second column must be IID. The order of individuals must be the same as val_pheno. Covariates that are factors (i.e. assessment center) must be encoded as such. The covariates are not penalized in the endoPRS model.
#' @param pheno_gwas A data frame containing the results of the GWAS run on the phenotype. It can be different than the one used in fit_endoPRS_single_iter as it can be run on the combined training and validation set. It must contain columns corresponding to chromosome, SNP, allele1, allele2, and P-value.
#' @param endo_gwas A data frame containing the results of the GWAS run on the endophenotype. It can be different than the one used in fit_endoPRS_single_iter as it can be run on the combined training and validation set. It must contain columns corresponding to chromosome, SNP, allele1, allele2, and P-value.
#' @param filter_hapmap An optional logical character. It corresponds to whether to only run the model on variants in the hapmap3 set.
#' @param hapmap An optional data frame with hapmap variants. This must be provided if filter_hapmap is set to TRUE. The hapmap3 data frame should contain columns CHR, SNP, A1, A2.
#' @param type An optional character vector of "linear" or "logistic." It specifies which type of model to fit to the data. If not provided, it will be learned from the phenotype. If the phenotype consists of 0's and 1's penalized logistic regression will be used, otherwie penalized linear regression will be used.
#' @param threshes An optional character of p-value thresholds. This corresponds to the GWAS p-value thresholds used to determine association with the phenotype and endophenotype. The optimal p-value threshold is learned from tuning on the validation set. If not provided, the original values from the endoPRS manuscript will be used (0.01, 1e-4, and 1e-6).
#' @param NCORES An optional value corresponding to the number of cores to use. Otherwise, the number of cores will be learned using nb_cores().
#' @param save_folder The path to a directory that contains the results of fit_endoPRS_single_iter when applied to the validation set.
#' @param save_model An optional boolean value. Indicates whether to save the final model and beta values to disk. If set to true, the save_folder will be used.
#'
#' @return A list with two elements: the betas of the final model and the model itself
#' \itemize{
#' \item{beta: A data frame consisting of the SNPs included in the final model and their corresponding coefficients. This can be used to apply the PRS to an external data set using software such as PLINK.}
#' \item{model: The final fitted endoPRS model as an object of class big_sp_list. This can be used to apply the final model to a test set using the predict() function from the bigsnpr package.  }
#' }
#'
#' @export
refit_endoPRS = function(G, map, fam,
                       train_pheno, train_covar,
                       val_pheno, val_covar,
                       pheno_gwas, endo_gwas,
                       filter_hapmap = F, hapmap = NULL, type = NULL,
                       threshes = c(0.01, 1e-4, 1e-6), NCORES = NULL,
                       save_folder,
                       save_model = T){

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


  ## Data frame to save validation results
  val_results = data.frame()

  ## Load the results of the models applied to validation
  for(thresh in threshes){
      for(iteration in 1:nrow(grid)){
        save_file = paste0(save_folder, "/trainingonly_model_thresh", thresh, "_w2", grid$w2[iteration],"_w3", grid$w3[iteration],"_validationresults.csv")
        res = read.csv(save_file)
        val_results = rbind(val_results, res)
      }
  }

  ## Determine the best performing model
  best_performing_model = val_results[which.max(val_res),]
  print(paste("Best Peforming Model Corresponds to Thresh:", best_performing_model$thresh,
              "w2:", best_performing_model$w2, "w3:", best_performing_model$w3 ))


  #######################################################################################
  ###########     Part 5: Refit using combined training + validation     ################
  #######################################################################################



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
  if(!is.null(save_model)){
    save_file = paste0(save_folder, "/combined_trainval_thresh", thresh, "_w2", w2,"_w3", w3,".rds")
    save(model, file = save_file)
  }

  ## Extract beta coefficients from fitted model
  beta_info = obtain_beta(model, ind.col, map)

  ## Return final model and final beta
  return(list(beta = beta_info, model = model))

}
