#' Function to fit weighted lasso model of endoPRS
#'
#' This function its a weighted lasso model. It penalizes the SNPs associated with only the phenotype, only the
#' endophenotype, and those associated with both differently. The model is either a penalized linear or penalized
#' logistic model.
#'
#'
#' @param G An object of class FBM from the bigsnpr package that contains genotypes.
#' @param y.train Vector of responses, correspoinding to train_index
#' @param train_index A vector of row indices that are used for training.
#' @param ind.col The vector of columns (SNPs) to use for fitting the model.
#' @param penalty_table The output of create_penalty_table a data frame that contains the name of the SNP, the group (pheno only, endo only, or both), and the associated penalty to be used in the weighted lasso model.
#' @param covar.train Matrix of covariates included in the model. They are unpenalized. They should correspond to train_index.
#' @param NCORES Number of cores to use for fitting the model.
#' @param type A character vector of "linear" or "logistic." It specifies which type of model to fit to the data.
#'
#' @return An object of class big_sp_list (a list of of 5) that contains the fitted model for each of the cross-validation splits.
#'

## Function to fit endoPRS lasso model
fit_weighted_model = function(G, y.train, train_index, ind.col, penalty_table, covar.train, NCORES, type){

  ## Run Linear Model
  if(type == "linear"){
    model <- big_spLinReg(X = G,
                          y.train,
                          ind.train = train_index,
                          ind.col = ind.col,
                          pf.X = penalty_table$penalty,
                          alphas = 1,
                          power_scale = 1,
                          K = 5,
                          nlambda = 100,
                          n.abort = 10,
                          dfmax = 100000,
                          covar.train = covar.train,
                          pf.covar = rep(0,ncol(covar.train)),
                          ncores = NCORES,
                          warn = F)

  }

  ## Run Logistic Model
  if(type == "logistic"){
    model <- big_spLogReg(X = G,
                          y.train,
                          ind.train = train_index,
                          ind.col = ind.col,
                          pf.X = penalty_table$penalty,
                          alphas = 1,
                          power_scale = 1,
                          K = 5,
                          nlambda = 100,
                          n.abort = 10,
                          dfmax = 100000,
                          covar.train = covar.train,
                          pf.covar = rep(0,ncol(covar.train)),
                          ncores = NCORES,
                          warn = F)

  }

  ## Return model
  return(model)
}
