library(pbmcapply)
library(dplyr)
library(parallel)
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

RunInteraction <- function(G2G_Obj,SOFTWARE_DIR,OUT_DIR,Ref_Panel,tool = 'PLINK',n_PC = 3,n_pPC = 3,
                           n_cores = 20,covars_to_incl = c(),lineage = c(),by_lineage = F,
                           debug=F,var_type = 'both',pheno = c('tb_score','xray_score'),X_Chr = T){
  PLINK = '/home/zmxu/Software/plink'
  GCTA = '/home/zmxu/Software/gcta'
  
  OUT_DIR <- glue::glue("{OUT_DIR}/{tool}/")
  system(glue::glue("mkdir -p {OUT_DIR}"))
  if(length(covars_to_incl) > 0){
    OUT_DIR <- glue::glue("{OUT_DIR}/PC_{n_PC}_pPC_{n_pPC}_cov_{paste0(sapply(covars_to_incl,function(x) gsub(x=x,pattern='_',replacement='')),collapse = '-')}/")
  }else{
    OUT_DIR <- glue::glue("{OUT_DIR}/PC_{n_PC}_pPC_{n_pPC}/")
  }
  system(glue::glue("mkdir -p {OUT_DIR}"))
  saveRDS(list(),glue::glue("{OUT_DIR}/xrayscore_int_results.rds"))
  saveRDS(list(),glue::glue("{OUT_DIR}/tbscore_int_results.rds"))
  
  if (tool == 'PLINK'){
    if(length(lineage)==0){
      lineages_to_run <- unique(names(G2G_Obj$aa_matrix_filt))
    }else{
      lineages_to_run <- lineage
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
      if(cur_lineage == 'ALL' & !by_lineage){
        cur_aa_matrix_filt <- G2G_Obj$aa_matrix_filt[[cur_lineage]]
      }else if(cur_lineage == 'ALL' & by_lineage){
        #SNP x LINEAGE interaction
        cur_aa_matrix_filt <- as.matrix(one_hot(as.data.table(data.frame(LINEAGE = G2G_Obj$vir_pPCs$ALL$LINEAGE))))
      }else{
        cur_aa_matrix_filt <- G2G_Obj$aa_matrix_filt[[cur_lineage]]
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
      tb_score <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),G2G_Obj$tb_score[[cur_lineage]] %>% dplyr::select(IID,tb_score),by=c('IID'='IID')) %>% dplyr::relocate(FID,IID,tb_score)
      data.table::fwrite(tb_score,col.names = T,quote = F,row.names = F,sep = ' ',file = glue::glue("{OUT_DIR_Lineage}/tmp/tb_score.txt"),na = 'NA')
      
      xray_score <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),G2G_Obj$xray_score[[cur_lineage]] %>% dplyr::select(IID,xray_score),by=c('IID'='IID')) %>% dplyr::relocate(FID,IID,xray_score)
      data.table::fwrite(xray_score,col.names = T,quote = F,row.names = F,sep = ' ',file = glue::glue("{OUT_DIR_Lineage}/tmp/xray_score.txt"),na = 'NA')
      
      #Run GWAS with lineage as covariate
      if(by_lineage){
        cur_covar <- cbind(covars[,1:2],data.frame(LINEAGE = G2G_Obj$vir_pPCs$ALL$LINEAGE),covars[,-c(1,2)])
        data.table::fwrite(cur_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_DIR_Lineage}/tmp/lineage_covar.txt"))
        
        lapply(pheno,function(x) system(
          glue::glue(
            "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --covar-variance-standardize --linear no-x-sex --pheno {OUT_DIR_Lineage}/tmp/{x}.txt --covar {OUT_DIR_Lineage}/tmp/lineage_covar.txt --out {OUT_DIR_Lineage}Lineage_GWAS"
          )
        ))
        
        discrete_covar <- dplyr::select(cur_covar,c('FID','IID',intersect(colnames(cur_covar),c('Patient_Sex','HIV_Status','LINEAGE'))))
        num_covar <- dplyr::select(cur_covar,c('FID','IID',setdiff(colnames(cur_covar),colnames(discrete_covar))))
        data.table::fwrite(discrete_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_DIR_Lineage}/tmp/lineage_discrete_covar.txt"))
        data.table::fwrite(num_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_DIR_Lineage}/tmp/lineage_num_covar.txt"))
        system(
          glue::glue(
            "{GCTA} --bfile {G2G_Obj$host_path[cur_lineage]} --make-grm --out {G2G_Obj$host_path[cur_lineage]} --thread-num {n_cores}"
          )
        )
        
        lapply(pheno,function(x) system(
          glue::glue(
            "{GCTA} --mlma --bfile {G2G_Obj$host_path[cur_lineage]} --grm {G2G_Obj$host_path[cur_lineage]} --autosome --pheno {OUT_DIR_Lineage}/tmp/{x}.txt --covar {OUT_DIR_Lineage}/tmp/lineage_discrete_covar.txt --qcovar {OUT_DIR_Lineage}/tmp/lineage_num_covar.txt --thread-num {n_cores} --out {OUT_DIR_Lineage}Lineage_GWAS"
          )
        ))
        
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
            
            lapply(pheno,function(x) system(
              glue::glue(
                "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --covar-variance-standardize --linear no-x-sex interaction --parameters 1-{ncol(cur_covar)} --pheno {OUT_DIR_Lineage}/tmp/{x}.txt --covar {OUT_DIR_Lineage}/tmp/{cur_pathogen_variant}.txt --out {OUT_DIR_Lineage}{cur_pathogen_variant}"
              )
            ))
          }else{
            cur_pathogen_variant <- colnames(cur_aa_matrix_filt)[k]
            cur_covar <- cbind(covars[,1:2],cur_aa_matrix_filt[,k,drop=FALSE],covars[,-c(1,2)])
            data.table::fwrite(cur_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_DIR_Lineage}/tmp/{cur_pathogen_variant}.txt"))
            
            if(X_Chr){
              lapply(pheno, function(x) system(
                glue::glue(
                  "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --covar-variance-standardize --linear interaction --parameters 1-{ncol(cur_covar)} --pheno {OUT_DIR_Lineage}/tmp/{x}.txt --covar {OUT_DIR_Lineage}/tmp/{cur_pathogen_variant}.txt --out {OUT_DIR_Lineage}{cur_pathogen_variant}"
                )
              ))
              
            }else{
              lapply(pheno, function(x) system(
                glue::glue(
                  "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --covar-variance-standardize --linear no-x-sex interaction --parameters 1-{ncol(cur_covar)} --pheno {OUT_DIR_Lineage}/tmp/{x}.txt --covar {OUT_DIR_Lineage}/tmp/{cur_pathogen_variant}.txt --out {OUT_DIR_Lineage}{cur_pathogen_variant}"
                )
              ))
            }
            
            lapply(pheno,function(x) system(glue::glue("grep 'ADD' {OUT_DIR_Lineage}{cur_pathogen_variant}.{x}.glm.linear > {OUT_DIR_Lineage}{cur_pathogen_variant}.{x}.glm.linear.add")))
            lapply(pheno,function(x) system(glue::glue("{SOFTWARE_DIR}pigz --fast {OUT_DIR_Lineage}{cur_pathogen_variant}.{x}.glm.linear.add")))
            lapply(pheno,function(x) system(glue::glue("rm {OUT_DIR_Lineage}{cur_pathogen_variant}.{x}.glm.linear")))
          }
        }
        )
      }
      results_tbscore <- list(GetResults(OUT_DIR_Lineage,suffix = 'tb_score.glm.linear.add.gz',p_thresh = 5e-8,n_cores = n_cores,is_interaction = T,is_ordinal = F,tool = tool))
      names(results_tbscore) <- cur_lineage
      all_res_tbscore <- readRDS(glue::glue("{OUT_DIR}/tbscore_int_results.rds"))
      all_res_tbscore <- c(all_res_tbscore,results_tbscore)
      saveRDS(all_res_tbscore,glue::glue("{OUT_DIR}/tbscore_int_results.rds"))
      
      results_xrayscore <- list(GetResults(OUT_DIR_Lineage,suffix = 'xray_score.glm.linear.add.gz',p_thresh = 5e-8,n_cores = n_cores,is_interaction = T,is_ordinal = F,tool = tool))
      names(results_xrayscore) <- cur_lineage
      all_res_xrayscore <- readRDS(glue::glue("{OUT_DIR}/xrayscore_int_results.rds"))
      all_res_xrayscore <- c(all_res_xrayscore,results_xrayscore)
      saveRDS(all_res_xrayscore,glue::glue("{OUT_DIR}/xrayscore_int_results.rds"))
      
      
    }
  }else if (tool == 'ORDINAL' | tool == 'ORDINAL_GxE'){
    if(length(lineage)==0){
      lineages_to_run <- unique(names(G2G_Obj$both_IDs_to_keep))
    }else{
      lineages_to_run <- lineage
    }
    for(i in 1:length(lineages_to_run)){
      cur_lineage <- lineages_to_run[i]
      OUT_DIR_Lineage <- glue::glue("{OUT_DIR}/LINEAGE_{cur_lineage}/")
      system(glue::glue("mkdir -p {OUT_DIR_Lineage}"))
      system(glue::glue("mkdir -p {OUT_DIR_Lineage}/tmp"))
      
      cur_wd <- getwd()
      setwd(OUT_DIR_Lineage)
      
      
      #Get IDD and FID, ensure they're in correct order
      fam_file <- data.table::fread(glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"))
      colnames(fam_file)[1:2] <- c('FID','IID')
      if(!all(fam_file$IID == G2G_Obj$both_IDs_to_keep[[cur_lineage]]$FAM_ID)){
        stop('Sample order incorrect')
      }
      
      #Store AA matrix
      if(cur_lineage == 'ALL' & !by_lineage){
        cur_aa_matrix_filt <- G2G_Obj$aa_matrix_raw
      }else if(cur_lineage == 'ALL' & by_lineage){
        #SNP x LINEAGE interaction
        cur_aa_matrix_filt <- as.matrix(one_hot(as.data.table(data.frame(LINEAGE = G2G_Obj$vir_pPCs$ALL$LINEAGE))))
      }else{
        cur_aa_matrix_filt <- G2G_Obj$aa_matrix_filt[[cur_lineage]]
      }
      
      #Store PCs and covars
      host_PCs <- G2G_Obj$host_PCs[[cur_lineage]]
      host_PCs <- dplyr::select(host_PCs,paste0('PC',1:n_PC))
      
      if(n_pPC != 0){
        pPCs <- G2G_Obj$vir_pPCs[[cur_lineage]]
        pPCs <-dplyr::select(pPCs,paste0('PC',1:n_pPC))
        colnames(pPCs) <- paste0('pPC',1:n_pPC)
      }
      if(length(covars_to_incl) > 0){
        covars_num_discrete <- G2G_Obj$covars[[cur_lineage]]
        covars_num_discrete <- dplyr::select(covars_num_discrete,covars_to_incl)
        if(n_pPC != 0){
          covars <- cbind(fam_file[,c(1,2)],host_PCs,pPCs,covars_num_discrete)
        }else{
          covars <- cbind(fam_file[,c(1,2)],host_PCs,covars_num_discrete)
        }
      }else{
        if(n_pPC != 0){
          covars <- cbind(fam_file[,c(1,2)],host_PCs,pPCs)
        }else{
          covars <- cbind(fam_file[,c(1,2)],host_PCs)
        }
      }
      
      #Store Phenotype
      pheno <- G2G_Obj$tb_score[[cur_lineage]]
      
      outcome_name <- setdiff(colnames(pheno),c('FID','IID'))
      covar_names <- paste(setdiff(colnames(covars),c('FID','IID')),collapse = ' + ')
      
      #Split genotype by chromosome
      lapply(1:22,function(q) system(glue::glue("{PLINK}2 --bfile {G2G_Obj$host_path[cur_lineage]} --set-all-var-ids @:#[b37]\\$r_\\$a --chr {q} --make-bed --out {OUT_DIR_Lineage}/tmp/chr{q}")))
      
      for(k in 1:ncol(cur_aa_matrix_filt)){
        print(k)
        cur_pathogen_variant <- colnames(cur_aa_matrix_filt)[k]
        cur_pathogen_dosage <- data.frame(AA = cur_aa_matrix_filt[,k,drop=T])
        cur_pheno <- cbind(data.frame(famid=covars[,1],perid=covars[,2],faid=0,moid=0),cur_pathogen_dosage,covars[,-c(1,2)],pheno[,-c(1,2)])
        write.csv(cur_pheno,quote = F,row.names = F,file = glue::glue("{OUT_DIR_Lineage}/tmp/{cur_pathogen_variant}.txt"))
        
        tryCatch({
          if(tool == 'ORDINAL'){
            mclapply(1:22,function(q){
              system(glue::glue("nice -n 19 julia -e 'using OrdinalGWAS; const datadir = \"{OUT_DIR_Lineage}/tmp/\"; ordinalgwas(@formula({outcome_name} ~ AA + {covar_names}), datadir * \"{cur_pathogen_variant}.txt\", datadir * \"chr{q}\",pvalfile=\"chr{q}\",testformula=@formula({outcome_name} ~ snp + snp & AA))'"))
            },mc.cores = n_cores,mc.silent = T,mc.preschedule = F)
            
          }else if(tool == 'ORDINAL_GxE'){
            mclapply(1:22,function(q){
              system(glue::glue("nice -n 19 julia -e 'using OrdinalGWAS; const datadir = \"{OUT_DIR_Lineage}/tmp/\"; ordinalgwas(@formula({outcome_name} ~ AA + {covar_names}), datadir * \"{cur_pathogen_variant}.txt\", datadir * \"chr{q}\",pvalfile=\"chr{q}\",analysistype = \"gxe\",e = :AA,test = :score)'"))
            },mc.cores = n_cores,mc.silent = T,mc.preschedule = F)
          }
          system(glue::glue("head -n 1 {OUT_DIR_Lineage}chr1 > {OUT_DIR_Lineage}{cur_pathogen_variant}.pval"))
          trash <- sapply(1:22,function(q) system(glue::glue("cat {OUT_DIR_Lineage}chr{q} | grep -v 'chr' >> {OUT_DIR_Lineage}{cur_pathogen_variant}.pval")))
          system(glue::glue("mv {OUT_DIR_Lineage}ordinalgwas.null.txt {OUT_DIR_Lineage}{cur_pathogen_variant}.null"))
          system(glue::glue("rm {OUT_DIR_Lineage}chr*"))
          system(glue::glue("pigz --fast {OUT_DIR_Lineage}{cur_pathogen_variant}.pval"))
        })
      }
    }
  }
}

SOFTWARE_DIR <- '/home/zmxu/Software/'
args <- commandArgs(trailingOnly = TRUE) 
G2G_Obj <- readRDS(args[[1]])
OUT_DIR <- args[[2]]
Ref_Panel <- args[[3]]
tool <- args[[4]]
n_PC <- as.numeric(args[[5]])
n_pPC <- as.numeric(args[[6]])
covars_to_incl <- args[[7]]
if(covars_to_incl == 'NA'){
  covars_to_incl <- c()
}else{
  covars_to_incl <- strsplit(covars_to_incl,split = ',')[[1]]
}
X_Chr <- as.logical(args[[8]])
n_cores <- as.numeric(args[[9]])

RunInteraction(G2G_Obj,SOFTWARE_DIR,OUT_DIR,
               Ref_Panel = Ref_Panel,tool = tool,
               n_cores = n_cores,n_PC = n_PC,n_pPC = n_pPC,
               covars_to_incl = covars_to_incl,
               debug = F,by_lineage = F,var_type = 'both',pheno = c('tb_score','xray_score'),X_Chr = X_Chr)