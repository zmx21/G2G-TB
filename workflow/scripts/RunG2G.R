library(pbmcapply)
library(dplyr)
library(data.table)
library(glue)
library(stringr)

GetResults <- function(OUT_DIR,suffix = 'glm.logistic.hybrid',p_thresh=5e-8,n_cores=5,is_interaction = F,is_ordinal = F,tool = 'PLINK'){
  all_files <- dir(OUT_DIR)
  result_files <- all_files[sapply(all_files,function(x) grepl(pattern = suffix,x=x))]
  if(!is_interaction){
    if((tool == 'PLINK' | tool == 'PLINK-FIRTH') & !grepl(x=suffix,pattern = 'gz')){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0("awk '{ if (NR == 1 || $13 <= ",p_thresh,") {print} }' ",OUT_DIR,x)),mc.cores = n_cores)
    }else if ((tool == 'PLINK' | tool == 'PLINK-FIRTH') & grepl(x=suffix,pattern = 'gz')){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0('zcat ',OUT_DIR,x," | awk '{ if (NR == 1 || $13 <= ",p_thresh,") {print} }' ")),mc.cores = n_cores)
    }else if(tool == 'GMMAT-SCORE'){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0('zcat ',OUT_DIR,x," | awk '{ if (NR == 1 || $11 <= ",p_thresh,") {print} }' ")),mc.cores = n_cores)
    }
  }else{
    if(is_ordinal){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0("zcat ",OUT_DIR,x," | awk -F ',' '{ if ($6 <= ",p_thresh,") {print} }' ")),mc.cores = n_cores)
    }else{
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0("zcat ",OUT_DIR,x," | grep 'ADDx' | awk '{ if ($12 <= ",p_thresh,") {print} }' ")),mc.cores = n_cores)
    }
  }
  names(results) <- result_files
  return(results)
}

RunG2G <- function(G2G_Obj,SOFTWARE_DIR,OUT_DIR,tool = 'PLINK',n_PC = 5,n_pPC = 6,n_cores = 20,covars_to_incl = c(),
                   lineage = c(),stratified = T,debug=F,chr = c(),test_snp = c(),model = NA,var_type = 'both'){
  system(glue::glue("mkdir -p {OUT_DIR}"))
  OUT_DIR <- glue::glue("{OUT_DIR}/{tool}/")
  system(glue::glue("mkdir -p {OUT_DIR}"))
  OUT_DIR <- glue::glue("{OUT_DIR}/PC_{n_PC}_pPC_{n_pPC}/")
  system(glue::glue("mkdir -p {OUT_DIR}"))
  # if(length(covars_to_incl) > 0){
  #   OUT_DIR <- glue::glue("{OUT_DIR}/PC_{n_PC}_pPC_{n_pPC}_cov_{paste0(sapply(covars_to_incl,function(x) gsub(x=x,pattern='_',replacement='')),collapse = '-')}/")
  # }else{
  #   OUT_DIR <- glue::glue("{OUT_DIR}/PC_{n_PC}_pPC_{n_pPC}/")
  # }
  OUT_DIR <- glue::glue("{OUT_DIR}/Stratified_{str_to_title(as.character(stratified))}/")
  system(glue::glue("mkdir -p {OUT_DIR}"))
  
  #Initialize results object
  saveRDS(list(),glue::glue("{OUT_DIR}/G2G_results.rds"))
  
  if (tool == 'PLINK' | tool == 'PLINK-FIRTH' | tool == 'HLA-PLINK' | tool == 'SAIGE' | tool == 'GMMAT-SCORE' | tool == 'GMMAT-WALD'){
    if(length(lineage)==0 & stratified){
      lineages_to_run <- unique(names(G2G_Obj$aa_matrix_filt))
    }else if(length(lineage) > 0 & stratified){
      lineages_to_run <- lineage
    }else if(!stratified){
      lineages_to_run <- 'ALL'
    }
    for(i in 1:length(lineages_to_run)){
      cur_lineage <- lineages_to_run[i]
      OUT_PATH_Lineage <- glue::glue("{OUT_DIR}/LINEAGE_{cur_lineage}/")
      system(glue::glue("mkdir -p {OUT_PATH_Lineage}"))
      system(glue::glue("mkdir -p {OUT_PATH_Lineage}/tmp"))
      
      #Get IDD and FID, ensure they're in correct order
      fam_file <- data.table::fread(glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"))
      colnames(fam_file)[1:2] <- c('FID','IID')
      if(!all(fam_file$IID == G2G_Obj$both_IDs_to_keep[[cur_lineage]]$FAM_ID)){
        stop('Sample order incorrect')
      }
      
      #Extract AA matrix for current lineage
      if(stratified){
        cur_aa_matrix_filt <- G2G_Obj$aa_matrix_filt[[cur_lineage]]
      }else{
        cur_aa_matrix_filt <- G2G_Obj$aa_matrix_full
      }

      #Extract/convert dosage
      #If burden, > 1 (more than one variant in gene) also as present
      #Var_Type = Both, treat homo and hetero calls as present
      #Var_Type = Homo, treat homo calls as present, hetero calls as missing (encoded as NA)
      if(grepl(pattern = 'Burden_True',x=OUT_DIR)){
        cur_aa_matrix_filt[cur_aa_matrix_filt > 1] <- 1
      }
      else if(var_type == 'both' | var_type == 'homo'){
        cur_aa_matrix_filt[cur_aa_matrix_filt==2] <- 1
      }else{
        stop('Invalid Var Type')
      }
      
      #Append IDs to matrix
      cur_aa_Matrix_mat <- cbind(fam_file[,c(1,2)],as.matrix(cur_aa_matrix_filt))
      colnames(cur_aa_Matrix_mat) <- c('FID','IID',colnames(cur_aa_matrix_filt))
      
      data.table::fwrite(cur_aa_Matrix_mat,col.names = T,row.names = F,sep = ' ',file = glue::glue("{OUT_PATH_Lineage}/tmp/AA_outcome.txt"),na = 'NA',quote = F)
      
      #Store PCs and covars
      host_PCs <- G2G_Obj$host_PCs[[cur_lineage]]
      host_PCs <- dplyr::select(host_PCs,'IID',paste0('PC',1:n_PC))
      
      if(n_pPC != 0){
        pPCs <- G2G_Obj$vir_pPCs[[cur_lineage]] %>% dplyr::inner_join(G2G_Obj$both_IDs_to_keep[[cur_lineage]]) %>% dplyr::rename(IID=FAM_ID)
        pPCs <-dplyr::select(pPCs,IID,paste0('PC',1:n_pPC))
        colnames(pPCs) <- c('IID',paste0('pPC',1:n_pPC))
      }
      cur_covars_to_incl <- covars_to_incl
      #Don't include lineage if only running in single lineage.
      if((nrow(str_locate_all(pattern = 'L',cur_lineage)[[1]]) == 1) & cur_lineage != 'ALL'){
        cur_covars_to_incl <- setdiff(covars_to_incl,'LINEAGE')
      }
      
      if(length(cur_covars_to_incl) > 0){
        covars_num_discrete <- G2G_Obj$covars[[cur_lineage]]
        #Add lineage as covariate if more than one lineage in current group
        if(stratified){
          if('LINEAGE' %in% cur_covars_to_incl & (nrow(str_locate_all(pattern = 'L',cur_lineage)[[1]]) > 1 | cur_lineage == 'ALL')){
            lineage_df <- G2G_Obj$both_IDs_to_keep[[cur_lineage]] %>% dplyr::select(IID=FAM_ID,LINEAGE)
            lineage_df$LINEAGE <- as.factor(lineage_df$LINEAGE)
            lineage_df <- one_hot(as.data.table(lineage_df),cols = 'LINEAGE')
            lineage_df <- lineage_df %>% dplyr::select(-colnames(lineage_df)[ncol(lineage_df)])
            
            covars_num_discrete <- covars_num_discrete %>% dplyr::left_join(lineage_df)
          }
        }

        covars_num_discrete <- dplyr::select(covars_num_discrete,'IID',contains(cur_covars_to_incl))
        covars_num_discrete[is.na(covars_num_discrete)] <- 'NONE'
        
        if(n_pPC != 0){
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(pPCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(covars_num_discrete,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }else{
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(covars_num_discrete,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }
      }else{
        if(n_pPC != 0){
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(pPCs,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }else{
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }
      }
      covars[covars==''] <- NA
      data.table::fwrite(covars,col.names = F,row.names = F,sep = '\t',file = glue::glue("{OUT_PATH_Lineage}/tmp/plink-covars.txt"),na = 'NA',quote = F)
      
      #Run association study for each AA variant in the current lineage
      AA_Matrix_No_ID <- cur_aa_Matrix_mat[,-c(1,2)]
      if(debug){
        end_index <- min(c(ncol(AA_Matrix_No_ID),5))
      }else{
        end_index <- ncol(AA_Matrix_No_ID)
      }
      
      pbmclapply(1:end_index,function(k){
        cur_pathogen_variant <- colnames(AA_Matrix_No_ID)[k]
        if(tool == 'PLINK'){
          if(is.na(model)){
              tryCatch(system(
                glue::glue(
                  "~/Software/plink2 --threads {min(c(n_cores,22))} --bfile {G2G_Obj$host_path[cur_lineage]} --glm cc-residualize hide-covar no-x-sex --covar-variance-standardize --1 --pheno {OUT_PATH_Lineage}/tmp/AA_outcome.txt --pheno-col-nums {k+2} --covar {OUT_PATH_Lineage}/tmp/plink-covars.txt --out {OUT_PATH_Lineage}{cur_pathogen_variant}"
                )
              ))
          }else{
            tryCatch(system(
              glue::glue(
                "~/Software/plink2 --threads {min(c(n_cores,22))} --bfile {G2G_Obj$host_path[cur_lineage]} --glm {model} cc-residualize hide-covar no-x-sex --covar-variance-standardize --1 --pheno {OUT_PATH_Lineage}/tmp/AA_outcome.txt --pheno-col-nums {k+2} --covar {OUT_PATH_Lineage}/tmp/plink-covars.txt --out {OUT_PATH_Lineage}{cur_pathogen_variant}.{model}"
              )
            ))
            
          }
          system(glue::glue('pigz --fast {OUT_PATH_Lineage}{cur_pathogen_variant}*.hybrid'))
          
        }else if(tool == 'HLA-PLINK'){
          Alleles_Path <- gsub(x=G2G_Obj$host_path[cur_lineage],pattern = 'Imputed',replacement = 'HLA_Alleles')
          AA_Path <- gsub(x=G2G_Obj$host_path[cur_lineage],pattern = 'Imputed',replacement = 'HLA_AA')
          system(glue::glue("mkdir -p {OUT_PATH_Lineage}/HLA_Allele/"))
          system(glue::glue("mkdir -p {OUT_PATH_Lineage}/HLA_AA/"))
          
          tryCatch(system(
            glue::glue(
              "~/Software/plink2 --threads {min(c(n_cores,22))} --bfile {Alleles_Path} --glm dominant hide-covar --1 --pheno {OUT_PATH_Lineage}/tmp/AA_outcome.txt --pheno-col-nums {k+2} --covar {OUT_PATH_Lineage}/tmp/plink-covars.txt --out {OUT_PATH_Lineage}/HLA_Allele/{cur_pathogen_variant}"
            )
          ))
          
          tryCatch(system(
            glue::glue(
              "~/Software/plink2 --threads {min(c(n_cores,22))} --bfile {AA_Path} --glm dominant hide-covar --1 --pheno {OUT_PATH_Lineage}/tmp/AA_outcome.txt --pheno-col-nums {k+2} --covar {OUT_PATH_Lineage}/tmp/plink-covars.txt --out {OUT_PATH_Lineage}/HLA_AA/{cur_pathogen_variant}"
            )
          ))
        }
        else if(tool == 'PLINK-FIRTH'){
          tryCatch(system(
            glue::glue(
              "~/Software/plink2 --threads {min(c(n_cores,22))} --bfile {G2G_Obj$host_path[cur_lineage]} --no-sex --glm firth hide-covar --1 --pheno {OUT_PATH_Lineage}/tmp/AA_outcome.txt --pheno-col-nums {k+2} --covar {OUT_PATH_Lineage}/tmp/plink-covars.txt --out {OUT_PATH_Lineage}{cur_pathogen_variant}"
            )
          ))
        }
      },mc.cores = max(floor(n_cores / 22),1))
      results <- list(GetResults(OUT_PATH_Lineage,suffix = 'glm.logistic.hybrid.gz',p_thresh=5e-8,n_cores=n_cores,is_interaction = F,is_ordinal = F,tool = tool))
      names(results) <- cur_lineage
      all_res <- readRDS(glue::glue("{OUT_DIR}/G2G_results.rds"))
      all_res <- c(all_res,results)
      saveRDS(all_res,glue::glue("{OUT_DIR}/G2G_results.rds"))
    }
  }
}

args <- commandArgs(trailingOnly = TRUE)
G2G_Obj <- readRDS(args[[1]])
OUT_DIR <- args[[2]]
tool <- args[[3]]
n_PC <- as.numeric(args[[4]])
n_pPC <- as.numeric(args[[5]])
covars_to_incl <- args[[6]]
if(covars_to_incl == 'NA'){
  covars_to_incl <- c()
}else{
  covars_to_incl <- strsplit(covars_to_incl,split = ',')[[1]]
}
stratified <- as.logical(args[[7]])
homo_only <- as.logical(args[[8]])
n_cores <- as.numeric(args[[9]])

# G2G_Obj <- readRDS('../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True/G2G_Obj.rds')
# OUT_DIR <- '../results/Burden_False_SIFT_False_Del_True_HomoOnly_True/'
# tool <- 'PLINK'
# n_PC <- 3
# n_pPC <- 0
# covars_to_incl <- c('patient_sex')
# if(covars_to_incl == 'NA'){
#   covars_to_incl <- c()
# }else{
#   covars_to_incl <- strsplit(covars_to_incl,split = ',')[[1]]
# }
# stratified <- F
# homo_only <- T
# n_cores <- 10

RunG2G(G2G_Obj,SOFTWARE_DIR,OUT_DIR,tool = tool,n_cores = n_cores,n_PC = n_PC,n_pPC = n_pPC,debug = F,covars_to_incl = covars_to_incl,
       model = NA,stratified = stratified,var_type = ifelse(homo_only,'homo','both'))