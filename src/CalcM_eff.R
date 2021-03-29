#Input: Rows are samples, columns are variables. 
CalcM_eff <- function(pheno_mat,method,prune_thresh = NA){
  #Correlation matrix of the phenotypes
  cor_mat <- cor(pheno_mat)
  #Eigenvalues of the correlation matrix
  eigen_val <- eigen(cor_mat)$values
  #Number of phenotypes
  M = ncol(pheno_mat)
  
  if(method == 'Nyholt'){
    #Effective number of tests 
    M_eff <- 1 + ((M - 1)*(1 - (var(eigen_val)/M)))
  }else if(method == 'Li_and_Ji'){
    M_eff <- 0
    #Sum across all eigenvalues
    for(i in 1:M){
      cur_eigen <- abs(eigen_val[i])
      M_eff <- M_eff + ifelse(cur_eigen >= 1,1,0) + (cur_eigen - floor(cur_eigen))
    }
  }else if(method == 'Prune'){
    if(!is.numeric(prune_thresh)){
      stop('Specific Prune Threshold')
    }
    locus_to_test <- colnames(pheno_mat)
    locus_to_keep <- c()
    locus_to_remove <- c()
    while(length(locus_to_test) > 0){
      cur_locus <- locus_to_test[1]
      cor_with_locus <- cor_mat[cur_locus,]^2
      locus_to_remove <- c(locus_to_remove,setdiff(names(cor_with_locus)[cor_with_locus >= prune_thresh],cur_locus))
      locus_to_keep <- c(locus_to_keep,cur_locus)
      locus_to_test <- setdiff(locus_to_test[-1],locus_to_remove)
    }
    M_eff <- length(locus_to_keep)
  }
  return(M_eff)
  
}
# CalcM_eff(G2G_Obj_AFGR_Tanz$aa_matrix_filt$L3,method = 'Prune',prune_thresh = 1)

