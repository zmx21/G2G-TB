library(dplyr)
library(parallel)

Map1KG <- function(bim_1KG,cur_lineage,scratch_dir){
  cur_lineage <- paste0('LINEAGE_',cur_lineage)
  
  if(!file.exists(glue::glue("{scratch_dir}{cur_lineage}/mapped_ids.txt"))){
    bim_1KG$V1 <- as.character(bim_1KG$V1)
    cur_host_bim <- data.table::fread(glue::glue("{scratch_dir}{cur_lineage}/TB_DAR_Imputed_G2G.bim"),header = F)
    jned_file <- dplyr::left_join(cur_host_bim,bim_1KG,by=c('V1'='V1','V4'='V4','V5'='V5','V6'='V6'))
    
    data.table::fwrite(jned_file %>% dplyr::select(ID_Host=V2.x,ID_1KG=V2.y),sep = '\t',
                       file = glue::glue("{scratch_dir}{cur_lineage}/mapped_ids.txt"))
  }
}
RunPASCAL <- function(res_obj,cur_lineage,scratch_dir,results_dir,PASCAL_path,filt = T,max_n_cores = 4,gene_scoring = 'max'){
  cur_lineage <- paste0('LINEAGE_',cur_lineage)
  cur_mapped_ids <- data.table::fread(glue::glue("{scratch_dir}{cur_lineage}/mapped_ids.txt"))
  
  ID_1KG <- cur_mapped_ids$ID_1KG
  ID_1KG[ID_1KG == ""] <- NA
  write(ID_1KG,glue::glue("{scratch_dir}{cur_lineage}/mapped_1KG_ids.txt"))
  remove(cur_mapped_ids)
  gc()
  
  if(!filt){
    cur_result_files <- names(res_obj[[1]])
  }else{
    if(length(res_obj[[1]][sapply(res_obj[[1]],nrow) > 0]) == 0){
      return(NA)
    }
    cur_result_files <- names(res_obj[[1]][sapply(res_obj[[1]],nrow) > 0])
  }
  n_cores <- min(c(max_n_cores,length(cur_result_files)))
  system(glue::glue("mkdir -p {results_dir}/{cur_lineage}/PASCAL/"))
  
  #lapply(cur_result_files,function(x){
  mclapply(cur_result_files,function(x){
    #Write summary stats
    if(grepl(pattern = 'glm.logistic.hybrid.gz',x=x)){
      system(glue::glue("zcat {results_dir}/{cur_lineage}/{x} | awk '{{if(NR > 1) print $13}}' > {results_dir}/{cur_lineage}/PASCAL/{gsub(x=x,pattern='.gz',replacement='.tmp')}"))
    }else{
      system(glue::glue("zcat {results_dir}/{cur_lineage}/{x} | grep 'ADDx' | awk '{{print $12}}' > {results_dir}/{cur_lineage}/PASCAL/{gsub(x=x,pattern='.gz',replacement='.tmp')}"))
    }
    system(glue::glue("paste {scratch_dir}{cur_lineage}/mapped_1KG_ids.txt {results_dir}/{cur_lineage}/PASCAL/{gsub(x=x,pattern='.gz',replacement='.tmp')} | grep -v 'NA' > {results_dir}/{cur_lineage}/PASCAL/{gsub(x=x,pattern='.gz',replacement='.pascal')}"))
    system(glue::glue("pigz --fast {results_dir}/{cur_lineage}/PASCAL/{gsub(x=x,pattern='.gz',replacement='.pascal')}"))
    system(glue::glue("rm {results_dir}/{cur_lineage}/PASCAL/{gsub(x=x,pattern='.gz',replacement='.tmp')}"))
    
    #Run PASCAL
    cur_wd = getwd()
    setwd(PASCAL_path)
    if(gene_scoring == 'max'){
      system(glue::glue("./Pascal --pval {results_dir}/{cur_lineage}/PASCAL/{gsub(x=x,pattern='.gz',replacement='.pascal.gz')} --outdir={results_dir}/{cur_lineage}/PASCAL/ --customdir=./resources/1kg_AFR/ --custom=AFR --runpathway=on --genescoring=max"))
      
    }else if(gene_scoring == 'sum'){
      system(glue::glue("./Pascal --pval {results_dir}/{cur_lineage}/PASCAL/{gsub(x=x,pattern='.gz',replacement='.pascal.gz')} --outdir={results_dir}/{cur_lineage}/PASCAL/ --customdir=./resources/1kg_AFR/ --custom=AFR --runpathway=on --genescoring=sum"))
      
    }else{
      stop('gene_scoring misspecified')
    }

    setwd(cur_wd)
  #}) 
  },mc.cores = n_cores,mc.preschedule = F)
  
}

# args <- commandArgs(trailingOnly = T)
# res_obj <- readRDS(args[[1]])
# bim_1KG <- args[[2]]
# scratch_dir <- args[[3]]
# results_dir <- args[[4]]

res_obj <- readRDS('/home/zmxu/G2G_TB/results/Burden_False_SIFT_False_Del_True_HomoOnly_True/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_0/tbscore_int_results.rds')
lineages <- names(res_obj)
bim_1KG <- data.table::fread('/home/zmxu/Software/PASCAL/resources/1kg_AFR/misc/AFR.bim',header = F)
scratch_dir <- '/home/zmxu/G2G_TB/scratch/Burden_False_SIFT_False_Del_True_HomoOnly_True/AFGR_Tanz_ChrX/'
results_dir <- '/home/zmxu/G2G_TB/results/Burden_False_SIFT_False_Del_True_HomoOnly_True/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_0/'
PASCAL_path <- '/home/zmxu/Software/PASCAL/'
# lapply(lineages,function(x) Map1KG(bim_1KG,x,scratch_dir))
mclapply(lineages,function(x) RunPASCAL(res_obj[x],x,scratch_dir,results_dir,PASCAL_path,filt = T),mc.cores = length(res_obj),mc.preschedule = F)


