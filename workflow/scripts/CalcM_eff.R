#Input: Rows are samples, columns are variables. 
CalcM_eff <- function(pheno_mat,method,prune_thresh = NA){
  #Correlation matrix of the phenotypes
  cor_mat <- cor(pheno_mat,use = 'complete.obs')
  
  if(method =='Nyholt' | method == 'Li_and_Ji' ){
    #Eigenvalues of the correlation matrix
    eigen_val <- eigen(cor_mat)$values
    #Number of phenotypes
    M = ncol(pheno_mat)
  }

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
#Group perfectly linked variants together
ConstructVariantSets <- function(pheno_mat){
  #Correlation matrix of the phenotypes
  cor_mat <- cor(pheno_mat,use = 'complete.obs')
  
  locus_to_test <- colnames(pheno_mat)
  variant_sets <- list()
  counter = 1
  while(length(locus_to_test) > 0){
    cur_locus <- locus_to_test[1]
    cor_with_locus <- cor_mat[cur_locus,]^2
    linked_variants <- names(cor_with_locus)[cor_with_locus == 1]
    variant_sets[[counter]] <- linked_variants
    locus_to_test <- setdiff(locus_to_test[-1],linked_variants)
    counter <- counter + 1
  }
  return(variant_sets)
}


G2G_Obj <- readRDS('../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds')
M_Eff = sapply(G2G_Obj$aa_matrix_filt,function(x) CalcM_eff(x,method = 'Prune',prune_thresh = 1))
M_Eff_full = CalcM_eff(G2G_Obj$aa_matrix_full,method = 'Prune',prune_thresh = 1)

Variant_Sets <- ConstructVariantSets(G2G_Obj$aa_matrix_full)


