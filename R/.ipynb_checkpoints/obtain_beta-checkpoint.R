#' Function to obtain beta from fitted lasso model.
#'
#' This function takes as input the fitted model from fit_weighted_model. It extracts the coefficients of the model
#' estimated in all five folds and averages them. Only non zero coefficients are returned.
#'
#' @param model a big_spLinReg/big_spLogReg model created by the fit_weighted_model function.
#' @param ind.col a vector containing the indices of SNPs from G (genotype matrix) included as input in the fit_weighted_model.
#' @param map data frame containing information about SNPs in G (genotype matrix).
#'
#' @return A data frame consisting of the SNPs included in the final model and their corresponding coefficients.
#'
#'
## Function to obtain beta from fitted model
obtain_beta = function(model, ind.col, map){

  ## Empty beta to hold covariates
  beta = 0

  ## Loop over each of the fitted folds in CMSA
  for(i in 1:length(model[[1]]) ){
    beta = beta + model[[1]][[i]]$beta[1:length(ind.col)]
  }
  ## Divide by number of folds
  beta = beta/length(model[[1]])

  ## Add final beta to map data frame and remove the entries that are zero
  beta_info = map[ind.col,]
  beta_info$beta = beta
  beta_info = beta_info[beta_info$beta != 0 ,]

  return(beta_info)
}
