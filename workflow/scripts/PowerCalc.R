library(genpwr)
library(ggplot2)
source('./CalcM_eff.R')
GetMAC <- function(AA_Matrix_filt){
  MAC <- apply(AA_Matrix_filt,2,function(x) {
    no_na <- x[!is.na(x)]
    if(length(no_na) == 0){
      return(0)
    }
    no_na_binary <- ifelse(no_na == 0,0,1) #Set to binary matrix (0 for absent, 1 more present - either homo or hetero call or burden)
    return(ifelse(length(unique(no_na_binary)) != 1,min(table(no_na_binary)),length(no_na_binary) - min(table(no_na_binary))))
  })
  return(MAC)
}

Var_Tbl <- readRDS('../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True/Mtb_Var_Tbl.rds')
Lineage <- data.table::fread('../../data/pheno/metadata_Sinergia_final_dataset_human_bac_genome_available_QCed.txt') %>% dplyr::select(G_NUMBER,Lineage)
Var_Tbl_MAC <- GetMAC(Var_Tbl)
Var_Tbl_MAC_By_Lineage <- sapply(1:ncol(Var_Tbl),function(x) max(apply(table(Var_Tbl[,x],Lineage$Lineage),2,min)))

Case_CTL_Ratio <- c(0.005,0.01,0.015,0.02,0.025)
MAC_Thresh <- round(c(0.005,0.01,0.015,0.02,0.025) * nrow(Var_Tbl))
Alpha <- 5e-8 / sapply(MAC_Thresh,function(x) CalcM_eff(Var_Tbl[,Var_Tbl_MAC > x & Var_Tbl_MAC_By_Lineage > x],method = 'Prune',prune_thresh = 1))

pw <- lapply(1:length(MAC_Thresh),function(i) genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=990, Case.Rate=Case_CTL_Ratio[i], k=NULL,
                  MAF=seq(0.05, 0.25, 0.01), OR=c(2,5,15,25),Alpha=Alpha[i],
                  True.Model=c( "Additive"), 
                  Test.Model=c("Additive")))
pw <- lapply(pw,function(x) {
  colnames(x)[9] <- 'Power'
  return(x)
})
pw <- do.call(rbind,pw)

p <- ggplot2::ggplot(aes(x=MAF,y=Power,color = factor(Case.Rate)),data = pw) + geom_line() + facet_grid(~OR,labeller = label_both)

