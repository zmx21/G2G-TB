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

RunInteraction <- function(G2G_Obj,OUT_DIR,tool = 'PLINK',n_PC = 3,n_pPC = 0,n_cores = 20,covars_to_incl = c(),lineage = c(),by_lineage = F,stratified = T,debug=F,model = NA,var_type = 'both',pheno = 'TB_score'){
  OUT_DIR <- glue::glue("{OUT_DIR}/{tool}/")
  system(glue::glue("mkdir -p {OUT_DIR}"))
  OUT_DIR <- glue::glue("{OUT_DIR}/PC_{n_PC}_pPC_{n_pPC}/")
  system(glue::glue("mkdir -p {OUT_DIR}"))
  OUT_DIR <- glue::glue("{OUT_DIR}/Stratified_{str_to_title(as.character(stratified))}/")
  system(glue::glue("mkdir -p {OUT_DIR}"))
  
  saveRDS(list(),glue::glue("{OUT_DIR}/{pheno}_int_results.rds"))

  if(length(lineage)==0 & stratified){
    lineages_to_run <- unique(names(G2G_Obj$aa_matrix_filt))
  }else if(length(lineage) > 0 & stratified){
    lineages_to_run <- lineage
  }else if(!stratified){
    lineages_to_run <- 'ALL'
  }
  
  for(i in 1:length(lineages_to_run)){
    cur_lineage <- lineages_to_run[i]
    OUT_DIR_Lineage <- glue::glue("{OUT_DIR}/LINEAGE_{cur_lineage}/")
    system(glue::glue("mkdir -p {OUT_DIR_Lineage}"))
    system(glue::glue("mkdir -p {OUT_DIR_Lineage}/tmp"))
    
    #Get IDD and FID, ensure they're in correct order
    fam_file <- data.table::fread(glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"))
    colnames(fam_file)[1:2] <- c('FID','IID')
    if(!all(fam_file$IID == G2G_Obj$both_IDs_to_keep[[cur_lineage]]$FAM_ID)){
      stop('Sample order incorrect')
    }
    
    #Store AA matrix
    if(!by_lineage & stratified){
      cur_aa_matrix_filt <- G2G_Obj$aa_matrix_filt[[cur_lineage]]
    }else if(!by_lineage & !stratified){
      cur_aa_matrix_filt <- G2G_Obj$aa_matrix_full
    }
    else if(by_lineage){
      #SNP x LINEAGE interaction
      cur_aa_matrix_filt <- as.matrix(one_hot(as.data.table(data.frame(LINEAGE = G2G_Obj$vir_pPCs$ALL$LINEAGE))))
    }
    
    #Extract/convert dosage
    #Var_Type = Both, treat homo and hetero calls as present
    #Var_Type = Homo, treat homo calls as present, hetero calls as absent
    if(var_type == 'both'){
      cur_aa_matrix_filt[cur_aa_matrix_filt==2] <- 1
    }else if (var_type == 'homo'){
      cur_aa_matrix_filt[cur_aa_matrix_filt==1] <- 0
      cur_aa_matrix_filt[cur_aa_matrix_filt==2] <- 1
    }else{
      stop('Invalid Var Type')
    }
    
    #Store PCs and covars
    host_PCs <- G2G_Obj$host_PCs[[cur_lineage]]
    host_PCs <- dplyr::select(host_PCs,'IID',paste0('PC',1:n_PC))
    
    if(n_pPC != 0){
      pPCs <- G2G_Obj$vir_pPCs[[cur_lineage]]
      pPCs <-dplyr::select(pPCs,IID,paste0('PC',1:n_pPC))
      colnames(pPCs) <- c('IID',paste0('pPC',1:n_pPC))
    }
    if(length(covars_to_incl) > 0){
      covars_num_discrete <- G2G_Obj$covars[[cur_lineage]]
      covars_num_discrete <- dplyr::select(covars_num_discrete,'IID',covars_to_incl)
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
    #Change categorial binary covariates to numeric (0 and 1)
    cat_covars <- which(sapply(3:ncol(covars),function(q) any(is.character(as.vector(t(covars[,..q]))))))
    if(length(cat_covars) > 0){
      for(col in (cat_covars+2)){
        covars[,col] <- as.integer(as.factor(t(covars[,..col]))) - 1
      }
    }
    
    #Write out Phenotype
    if(pheno == 'TB_score'){
      tb_score <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),G2G_Obj$tb_score[[cur_lineage]] %>% dplyr::select(IID,TB_score),by=c('IID'='IID')) %>% dplyr::relocate(FID,IID,TB_score)
      data.table::fwrite(tb_score,col.names = T,quote = F,row.names = F,sep = ' ',file = glue::glue("{OUT_DIR_Lineage}/tmp/TB_score.txt"),na = 'NA')
    }else if(pheno == 'Xray_score'){
      xray_score <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),G2G_Obj$xray_score[[cur_lineage]] %>% dplyr::select(IID,Xray_score),by=c('IID'='IID')) %>% dplyr::relocate(FID,IID,Xray_score)
      data.table::fwrite(xray_score,col.names = T,quote = F,row.names = F,sep = ' ',file = glue::glue("{OUT_DIR_Lineage}/tmp/Xray_score.txt"),na = 'NA')
      
    }else{
      stop('Invalid Pheno')
    }


    #Run GWAS with lineage as covariate
    if(by_lineage){
      cur_covar <- cbind(covars[,1:2],data.frame(LINEAGE = G2G_Obj$vir_pPCs$ALL$LINEAGE),covars[,-c(1,2)])
      data.table::fwrite(cur_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_DIR_Lineage}/tmp/lineage_covar.txt"))
      
      system(
        glue::glue(
          "plink2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --covar-variance-standardize --linear no-x-sex --pheno {OUT_DIR_Lineage}/tmp/{pheno}.txt --covar {OUT_DIR_Lineage}/tmp/lineage_covar.txt --out {OUT_DIR_Lineage}Lineage_GWAS"
        )
      )
      
      discrete_covar <- dplyr::select(cur_covar,c('FID','IID',intersect(colnames(cur_covar),c('Patient_Sex','HIV_Status','LINEAGE'))))
      num_covar <- dplyr::select(cur_covar,c('FID','IID',setdiff(colnames(cur_covar),colnames(discrete_covar))))
      data.table::fwrite(discrete_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_DIR_Lineage}/tmp/lineage_discrete_covar.txt"))
      data.table::fwrite(num_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_DIR_Lineage}/tmp/lineage_num_covar.txt"))
      system(
        glue::glue(
          "gcta64 --bfile {G2G_Obj$host_path[cur_lineage]} --make-grm --out {G2G_Obj$host_path[cur_lineage]} --thread-num {n_cores}"
        )
      )
      
      system(
        glue::glue(
          "gcta64 --mlma --bfile {G2G_Obj$host_path[cur_lineage]} --grm {G2G_Obj$host_path[cur_lineage]} --autosome --pheno {OUT_DIR_Lineage}/tmp/{pheno}.txt --covar {OUT_DIR_Lineage}/tmp/lineage_discrete_covar.txt --qcovar {OUT_DIR_Lineage}/tmp/lineage_num_covar.txt --thread-num {n_cores} --out {OUT_DIR_Lineage}Lineage_GWAS"
        )
      )
    }
    if(debug){
      end_index <- min(c(ncol(cur_aa_matrix_filt),2))
    }else{
      end_index <- ncol(cur_aa_matrix_filt)
    }
    
    for(k in 1:end_index){
      tryCatch({
        if(by_lineage){
          cur_pathogen_variant <- colnames(cur_aa_matrix_filt)[k]
          cur_covar <- cbind(covars[,1:2],cur_aa_matrix_filt[,k,drop=FALSE],cur_aa_matrix_filt[,-k][,-1],covars[,-c(1,2)])
          data.table::fwrite(cur_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_DIR_Lineage}/tmp/{cur_pathogen_variant}.txt"))
          
          system(
            glue::glue(
              "plink2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --covar-variance-standardize --linear no-x-sex interaction --parameters 1-{ncol(cur_covar)} --pheno {OUT_DIR_Lineage}/tmp/{pheno}.txt --covar {OUT_DIR_Lineage}/tmp/{cur_pathogen_variant}.txt --out {OUT_DIR_Lineage}{cur_pathogen_variant}"
            )
          )
        }else{
          cur_pathogen_variant <- colnames(cur_aa_matrix_filt)[k]
          cur_covar <- cbind(covars[,1:2],cur_aa_matrix_filt[,k,drop=FALSE],covars[,-c(1,2)])
          data.table::fwrite(cur_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_DIR_Lineage}/tmp/{cur_pathogen_variant}.txt"))
          
          system(
            glue::glue(
              "plink2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --covar-variance-standardize --linear no-x-sex interaction --parameters 1-{ncol(cur_covar)} --pheno {OUT_DIR_Lineage}/tmp/{pheno}.txt --covar {OUT_DIR_Lineage}/tmp/{cur_pathogen_variant}.txt --out {OUT_DIR_Lineage}{cur_pathogen_variant}"
            )
          )
          
          
          system(glue::glue("grep 'ADD' {OUT_DIR_Lineage}{cur_pathogen_variant}.{pheno}.glm.linear > {OUT_DIR_Lineage}{cur_pathogen_variant}.{pheno}.glm.linear.add"))
          system(glue::glue("pigz --fast {OUT_DIR_Lineage}{cur_pathogen_variant}.{pheno}.glm.linear.add"))
          system(glue::glue("rm {OUT_DIR_Lineage}{cur_pathogen_variant}.{pheno}.glm.linear"))
        }
      }
      )
    }
    results <- list(GetResults(OUT_DIR_Lineage,suffix = glue::glue('{pheno}.glm.linear.add.gz'),p_thresh = 5e-8,n_cores = n_cores,is_interaction = T,is_ordinal = F,tool = tool))
    names(results) <- cur_lineage
    all_res <- readRDS(glue::glue("{OUT_DIR}/{pheno}_int_results.rds"))
    all_res <- c(all_res,results)
    saveRDS(all_res,glue::glue("{OUT_DIR}/{pheno}_int_results.rds"))
    
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
pheno <- args[[9]]
n_cores <- as.numeric(args[[10]])

# G2G_Obj <- readRDS('../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True/G2G_Obj.rds')
# OUT_DIR <- '../results/Burden_False_SIFT_False_Del_False_HomoOnly_True/'
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
# pheno <- 'TB_score'

RunInteraction(G2G_Obj,OUT_DIR,
               tool = tool,
               n_cores = n_cores,n_PC = n_PC,n_pPC = n_pPC,
               covars_to_incl = covars_to_incl,
               debug = F,by_lineage = F,stratified=stratified,var_type = ifelse(homo_only,'homo','both'),pheno = pheno)