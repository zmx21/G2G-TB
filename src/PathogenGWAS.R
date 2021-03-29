library(MASS)
library(pbmcapply)
library(dplyr)
library(broom)
PathogenGWAS <- function(G2G_Obj){
  pgwas_result <- list()
  for(i in 1:length(G2G_Obj$aa_matrix_filt)){
    cur_mat <- G2G_Obj$aa_matrix_filt[[i]]
    metadata <- G2G_Obj$covars[[i+1]] %>% dplyr::select(-FID) 
    
    result <- pbmclapply(colnames(cur_mat),function(x){
      cur_variant <- cur_mat[,x]
      cur_data <- data.frame(Patho_SNP = cur_variant,
                             metadata,
                             tb_score = factor(G2G_Obj$tb_score[[i+1]]$tb_score[match(G2G_Obj$both_IDs_to_keep[[i+1]]$FAM_ID,G2G_Obj$tb_score[[i+1]]$IID)]))
      fit <- clm(glue::glue("tb_score ~ Patho_SNP + Age + Patient_Sex"),data = cur_data)
      fit_summary <- summary(fit)
      return(fit_summary$coefficients)
    },mc.cores = 3)    
    pgwas_result <- c(pgwas_result,list(result))
  }
  names(pgwas_result) <- names(G2G_Obj$aa_matrix_filt)
  return(pgwas_result)
}
G2G_Obj_AFGR_Tanz <- readRDS('~/G2G_TB/scratch/AFGR_Tanz/G2G_Obj_AFGR_Tanz.rds')
pgwas_result <- PathogenGWAS(G2G_Obj_AFGR_Tanz)

p_L1 <- sapply(pgwas_result$L1,function(x) x['Patho_SNP','Pr(>|z|)'])
p_L2 <- sapply(pgwas_result$L2,function(x) x['Patho_SNP','Pr(>|z|)'])
p_L3 <- sapply(pgwas_result$L3,function(x) x['Patho_SNP','Pr(>|z|)'])
p_L4 <- sapply(pgwas_result$L4,function(x) x['Patho_SNP','Pr(>|z|)'])
