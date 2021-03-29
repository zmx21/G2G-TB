library(ordinal)
library(dplyr)
library(broom)
library(ggplot2)
library(gtsummary)
library(data.table)
library(pbmcapply)
my_tidy <- function(x, exponentiate =  FALSE, conf.level = 0.95, ...) {
  dplyr::bind_cols(
    broom::tidy(x, exponentiate = exponentiate, conf.int = FALSE),
    # calculate the confidence intervals, and save them in a tibble
    stats::confint.default(x) %>%
      tibble::as_tibble() %>%
      rlang::set_names(c("conf.low", "conf.high"))  )
}

FitOrdinal <- function(G2G_Obj,host_snp,lineage = 'ALL'){
  host_genotype <- snpStats::read.plink(bed = G2G_Obj$host_path[[lineage]],select.snps = host_snp)
  host_dosage <- as(host_genotype$genotypes, Class = 'numeric')
  host_dosage_ordered <- host_dosage[match(G2G_Obj$both_IDs_to_keep[[lineage]]$FAM_ID,rownames(host_dosage)),]

  metadata <- data.frame(IID = names(host_dosage_ordered),Host_SNP = host_dosage_ordered,
                         tb_score = factor(G2G_Obj$tb_score[[lineage]]$tb_score[match(G2G_Obj$both_IDs_to_keep[[lineage]]$FAM_ID,G2G_Obj$tb_score[[lineage]]$IID)]))
  
  metadata <- dplyr::left_join(metadata,by=c('IID' = 'IID'),G2G_Obj$covars[[lineage]] %>% dplyr::select(-FID)) %>%
    dplyr::left_join(G2G_Obj$host_PCs[[lineage]],by=c('IID' = 'IID')) %>%
    dplyr::left_join(G2G_Obj$both_IDs_to_keep[[lineage]] %>% dplyr::select(IID=FAM_ID,LINEAGE),by=c('IID' = 'IID'))
  lineage_covar <- mltools::one_hot(as.data.table(data.frame(LINEAGE = metadata$LINEAGE)))
  metadata <- cbind(metadata,lineage_covar)
  fit <- clm("tb_score ~  Host_SNP + Patient_Sex + Age + PC1 + PC2 + PC3 + LINEAGE_L3 + LINEAGE_L2 + LINEAGE_L4 + HIV_Status",data = metadata)
  fit_summary <- summary(fit)
  return(fit_summary$coefficients['Host_SNP','Pr(>|z|)'])
}

FitOrdinalInt <- function(G2G_Obj,host_snp,AA_variant,lineage = NA,by_lineage = F,plot = F){
  host_genotype <- snpStats::read.plink(bed = G2G_Obj$host_path["ALL"],select.snps = host_snp)
  host_dosage <- as(host_genotype$genotypes, Class = 'numeric')

  if(by_lineage){
    host_dosage_ordered <- host_dosage[match(G2G_Obj$both_IDs_to_keep[["ALL"]]$FAM_ID,rownames(host_dosage)),]
    pathogen_dosage_ordered <- factor(G2G_Obj$vir_pPCs$ALL$LINEAGE == lineage)
    names(pathogen_dosage_ordered) <- G2G_Obj$both_IDs_to_keep[["ALL"]]$G_NUMBER
    metadata <- data.frame(IID = names(host_dosage_ordered),G_NUMBER = names(pathogen_dosage_ordered),Host_SNP = host_dosage_ordered,Patho_Variant = pathogen_dosage_ordered,
                           tb_score = factor(G2G_Obj$tb_score$ALL$tb_score[match(G2G_Obj$both_IDs_to_keep[["ALL"]]$FAM_ID,G2G_Obj$tb_score$ALL$IID)]))
    
  }else{
    if(!is.na(lineage)){
      host_dosage_ordered <- host_dosage[match(G2G_Obj$both_IDs_to_keep[[lineage]]$FAM_ID,rownames(host_dosage)),]
      pathogen_dosage <- G2G_Obj$aa_matrix_filt[[lineage]][,AA_variant]
      pathogen_dosage_ordered <- pathogen_dosage[match(G2G_Obj$both_IDs_to_keep[[lineage]]$G_NUMBER,names(pathogen_dosage))]
      metadata <- data.frame(IID = names(host_dosage_ordered),G_NUMBER = names(pathogen_dosage_ordered),Host_SNP = host_dosage_ordered,Patho_Variant = pathogen_dosage_ordered,
                             tb_score = factor(G2G_Obj$tb_score[[lineage]]$tb_score[match(G2G_Obj$both_IDs_to_keep[[lineage]]$FAM_ID,G2G_Obj$tb_score[[lineage]]$IID)]))
      
    }else{
      host_dosage_ordered <- host_dosage[match(G2G_Obj$both_IDs_to_keep[["ALL"]]$FAM_ID,rownames(host_dosage)),]
      pathogen_dosage <- G2G_Obj$aa_matrix_raw[,AA_variant]
      pathogen_dosage_ordered <- pathogen_dosage[match(G2G_Obj$both_IDs_to_keep[["ALL"]]$G_NUMBER,names(pathogen_dosage))]
      metadata <- data.frame(IID = names(host_dosage_ordered),G_NUMBER = names(pathogen_dosage_ordered),Host_SNP = host_dosage_ordered,Patho_Variant = pathogen_dosage_ordered,
                             tb_score = factor(G2G_Obj$tb_score$ALL$tb_score[match(G2G_Obj$both_IDs_to_keep[["ALL"]]$FAM_ID,G2G_Obj$tb_score$ALL$IID)]))
      
    }
  }
  
  metadata <- dplyr::left_join(metadata,by=c('IID' = 'IID'),G2G_Obj$covars$ALL %>% dplyr::select(-FID)) %>%
    dplyr::left_join(G2G_Obj$host_PCs$ALL,by=c('IID' = 'IID')) %>%
    dplyr::left_join(G2G_Obj$vir_pPCs$ALL %>% dplyr::select(G_NUMBER,pPC1 = PC1,pPC2 = PC2, pPC3 = PC3,LINEAGE),by=c('G_NUMBER' = 'G_NUMBER'))

  if(by_lineage){
    lineage_covar <- mltools::one_hot(as.data.table(data.frame(LINEAGE = metadata$LINEAGE)))
    metadata <- cbind(metadata,lineage_covar)
    fit <- clm(glue::glue("tb_score ~ LINEAGE_{lineage} * Host_SNP + Patient_Sex + Age + PC1 + PC2 + PC3 + {paste0('LINEAGE_',setdiff(c('L1','L2','L3','L4'),lineage))[-1]} + HIV_Status"),data = metadata)
    fit_summary <- summary(fit)
    p_val <- fit_summary$coefficients[paste0('LINEAGE_',lineage,':Host_SNP'),'Pr(>|z|)']
    if(plot){
      metadata$tb_score <- as.numeric(metadata$tb_score)
      formula_lm <- paste0("tb_score ~ Patient_Sex + Age + HIV_Status + PC1 + PC2 + PC3")
      fit_lm <- lm(formula_lm,data = metadata)
      metadata$corrected_tb <- rep(NA,nrow(metadata))
      corrected_tb <- residuals(fit_lm)
      metadata$corrected_tb[match(names(corrected_tb),rownames(metadata))] <- corrected_tb
      
      ggplot2::ggplot(data = metadata,aes(x = factor(Host_SNP),y=corrected_tb)) + geom_boxplot() + facet_grid(~factor(LINEAGE_L3),labeller = label_both) + xlab(host_snp) + ylab('Corrected TB Score') 
    }
  }else{
    fit <- clm('tb_score ~ Patho_Variant * Host_SNP  + Patient_Sex + Age + PC1 + PC2 + PC3',data=metadata)
    fit_summary <- summary(fit)
    p_val <- fit_summary$coefficients['Patho_Variant:Host_SNP','Pr(>|z|)']
  }
  return(p_val)
}
G2G_Obj_AFGR_Tanz <- readRDS('~/G2G_TB/scratch/AFGR_Tanz/G2G_Obj_AFGR_Tanz.rds')
for(lineage in c('L1','L2','L3','L4')){
  by_lineage <- data.table::fread(glue::glue("~/G2G_TB/interaction_results/AFGR_Tanz/PLINK/PC_3_pPC_0_cov_PatientSex-Age-HIVStatus/LINEAGE_ALL/LINEAGE_{lineage}.tb_score.glm.linear"))  %>%
    dplyr::filter(grepl(x=TEST,pattern = 'ADDx')) %>% dplyr::arrange(P)
  by_lineage_filt <- dplyr::filter(by_lineage, P < 1e-3)
  remove(by_lineage)
  if(nrow(by_lineage_filt) == 0){
    next
  }
  p_res <- pbmclapply(by_lineage_filt$ID,function(x) FitOrdinalInt(G2G_Obj_AFGR_Tanz,host_snp = x,lineage = lineage,by_lineage = T),mc.cores = 10)
  by_lineage_filt$P_Ordinal <- unlist(p_res)
  data.table::fwrite(by_lineage_filt,glue::glue("~/G2G_TB/interaction_results/AFGR_Tanz/PLINK/PC_3_pPC_0_cov_PatientSex-Age-HIVStatus/LINEAGE_ALL/LINEAGE_{lineage}.tb_score.ordinal"),sep = ' ',col.names = T,row.names = F,na = 'NA')
}

# L3_res <- readRDS('/home/zmxu/G2G_TB/interaction_results/AFGR_Tanz/PLINK/PC_3_pPC_0_cov_PatientSex-Age-HIVStatus/LINEAGE_L3/results.rds')
# L3_res_filt <- do.call(rbind,L3_res[sapply(L3_res,nrow) > 0])
# p_res <- pbmclapply(1:nrow(L3_res_filt),function(i) FitOrdinalInt(G2G_Obj_AFGR_Tanz,host_snp = L3_res_filt$V3[i],AA_variant =  gsub(pattern = 'ADDx',replacement = '',x = L3_res_filt$V7[i]),lineage = 'L3',by_lineage = F),mc.cores = 10)
# L3_res_filt$p_ordinal <- p_res
# data.table::fwrite(L3_res_filt,'/home/zmxu/G2G_TB/interaction_results/AFGR_Tanz/PLINK/PC_3_pPC_0_cov_PatientSex-Age-HIVStatus/LINEAGE_L3/ordinal.result.txt',col.names = F,row.names = F,sep = ' ',na = 'NA')



# L4_res <- readRDS('/home/zmxu/G2G_TB/interaction_results/AFGR_Tanz/PLINK/PC_3_pPC_0_cov_PatientSex-Age-HIVStatus/LINEAGE_L4/results.rds')
# L4_res_filt <- do.call(rbind,L4_res[sapply(L4_res,nrow) > 0])
# p_res <- pbmclapply(1:nrow(L4_res_filt),function(i) FitOrdinalInt(G2G_Obj_AFGR_Tanz,host_snp = L4_res_filt$V3[i],AA_variant =  gsub(pattern = 'ADDx',replacement = '',x = L4_res_filt$V7[i]),lineage = 'L4',by_lineage = F),mc.cores = 10)
# L4_res_filt$p_ordinal <- p_res
# data.table::fwrite(L4_res_filt,'/home/zmxu/G2G_TB/interaction_results/AFGR_Tanz/PLINK/PC_3_pPC_0_cov_PatientSex-Age-HIVStatus/LINEAGE_L4/ordinal.result.txt',col.names = F,row.names = F,sep = ' ',na = 'NA')

# FitOrdinal(G2G_Obj_AFGR_Tanz,host_snp = 'rs79001959')