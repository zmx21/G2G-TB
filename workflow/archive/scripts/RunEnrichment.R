library(pbmcapply)
ParseMtbNetwork <- function(tbl_path){
  #Genes of each cluster
  tbl <- data.table::fread(tbl_path)
  genes <- lapply(tbl$genes,function(x) unlist(strsplit(x=x,split = '\\|')))
  names(genes) <- tbl$bicluster
  
  #Regulators of each cluster 
  regulators <- lapply(tbl$regulator,function(x) strsplit(gsub(x=gsub(x=x,pattern = '-up-regulates-bicluster',replacement = ''),pattern = '-down-regulates-bicluster',replacement = ''),split = '\\|')[[1]])
  regulators <- lapply(regulators,function(x) {
    if(length(x) == 1){
      if(is.na(x)){
        return(NA)
      }
    }else{
     return(sapply(x,function(y) strsplit(x=y,split = '_')[[1]][1],USE.NAMES = F)) 
    }})
  names(regulators) <- tbl$bicluster
  #Info for each cluster
  tbl_info <- dplyr::select(tbl,bicluster,n_genes=nrow,score,functions)
  
  return(list(genes=genes,regulators=regulators,tbl_info=tbl_info))
}
GetMaxP <- function(results_file){
  res <- data.table::fread(cmd=glue::glue("zcat {results_file} | grep \'x\' | awk \'{{print $12}}\'"),header = F)
  min_p <- min(res$V1)
  names(min_p) <- strsplit(x=strsplit(results_file,split = '/')[[1]][length(strsplit(results_file,split = '/')[[1]])],split = ':')[[1]][1]
  return(min_p)
}

# mtb_network_tbl <- ParseMtbNetwork('~/G2G_TB/data/Mtb/MTB_Network/mtb_bicluster_attributes.txt')
# G2G_Res <- readRDS('~/G2G_TB/results/Burden_False_SIFT_False_Del_True_HomoOnly_True/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_0/G2G_results.rds')
# G2G_Res_Filt <- unlist(G2G_Res,recursive = F)
# G2G_Res_Filt <- G2G_Res_Filt[sapply(G2G_Res_Filt,nrow) > 0]
# 
# sig_genes <- sapply(names(G2G_Res_Filt),function(x) strsplit(strsplit(strsplit(x,split = '\\.')[[1]][2],split = ':')[[1]][1],split = '_')[[1]][2])

module_p <- list()
sig_in_module <- rep(NA,length(mtb_network_tbl$genes))
for(i in 1:length(mtb_network_tbl$genes)){
  sig_in_module[i] <- length(intersect(sig_genes,mtb_network_tbl$genes[[i]]))
  sig_not_in_module <- length(setdiff(sig_genes,mtb_network_tbl$genes[[i]]))
  not_sig_in_module <- length(mtb_network_tbl$genes[[i]]) - sig_in_module[i]
  not_sig_not_in_module <- length(unlist(mtb_network_tbl$genes)) - sig_in_module[i] - sig_not_in_module - not_sig_in_module
  module_p[[i]] <- fisher.test(matrix(c(sig_in_module[i],not_sig_in_module,sig_not_in_module,not_sig_not_in_module),ncol=2))
}



sapply(mtb_network_tbl$genes,function(x) length(intersect(x,sig_genes))) / sapply(mtb_network_tbl$genes,length)
# lineages <- paste0('~/G2G_TB/results/Burden_False_SIFT_False_Del_True_HomoOnly_True/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_0/',dir('~/G2G_TB/results/Burden_False_SIFT_False_Del_True_HomoOnly_True/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_0/')[grepl(dir('~/G2G_TB/results/Burden_False_SIFT_False_Del_True_HomoOnly_True/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_0/'),pattern = 'LINEAGE')])
# max_p <- lapply(lineages,function(x) unlist(pbmclapply(paste0(x,'/',dir(x)[grepl(x=dir(x),pattern = 'tb_score')]),function(y) GetMaxP(y),mc.cores = 50)))
# saveRDS(max_p,file = '~/G2G_TB/results/Burden_False_SIFT_False_Del_True_HomoOnly_True/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_0/max_p.rds')

