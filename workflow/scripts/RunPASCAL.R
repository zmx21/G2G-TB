library(dplyr)
library(parallel)

Map1KG <- function(bim_1KG,cur_lineage,scratch_dir){
  cur_lineage <- paste0('LINEAGE_',cur_lineage)
  
  if(!file.exists(glue::glue("{scratch_dir}{cur_lineage}/mapped_ids.txt"))){
    bim_1KG <- data.table::fread(bim_1KG,header = F)
    bim_1KG$V1 <- as.character(bim_1KG$V1)
    cur_host_bim <- data.table::fread(glue::glue("{scratch_dir}{cur_lineage}/TB_DAR_G2G.bim"),header = F)
    jned_file <- dplyr::left_join(cur_host_bim,bim_1KG,by=c('V1'='V1','V4'='V4','V5'='V5','V6'='V6'))
    
    data.table::fwrite(jned_file %>% dplyr::select(ID_Host=V2.x,ID_1KG=V2.y),sep = '\t',
                       file = glue::glue("{scratch_dir}{cur_lineage}/mapped_ids.txt"))
    
    cur_mapped_ids <- data.table::fread(glue::glue("{scratch_dir}{cur_lineage}/mapped_ids.txt"))
    
    ID_1KG <- cur_mapped_ids$ID_1KG
    ID_1KG[ID_1KG == ""] <- NA
    write(ID_1KG,glue::glue("{scratch_dir}{cur_lineage}/mapped_1KG_ids.txt"))
    remove(cur_mapped_ids)
    gc()
    
  }
}
RunPASCAL <- function(GWAS_file,variant,scratch_dir,results_dir,PASCAL_path,filt = T,gene_scoring = 'max',cur_lineage='ALL'){
  system(glue::glue("mkdir -p {results_dir}"))
  cur_lineage <- paste0('LINEAGE_',cur_lineage)

  #Write summary stats
  if(grepl(pattern = 'glm.logistic.hybrid.gz',x=GWAS_file)){
    system(glue::glue("zcat {GWAS_file} | awk '{{if(NR > 1) print $13}}' > {results_dir}/{variant}.tmp"))
  }else{
    system(glue::glue("zcat {GWAS_file} | grep 'ADDx' | awk '{{print $12}}' > {results_dir}/{variant}.tmp"))
  }
  system(glue::glue("paste {scratch_dir}{cur_lineage}/mapped_1KG_ids.txt {results_dir}/{variant}.tmp | grep -v 'NA' > {results_dir}/{variant}.pascal"))
  system(glue::glue("pigz --fast {results_dir}/{variant}.pascal"))
  system(glue::glue("rm {results_dir}/{variant}.tmp"))
  
  #Run PASCAL
  cur_wd = getwd()
  setwd(PASCAL_path)
  if(gene_scoring == 'max'){
    system(glue::glue("./Pascal --pval {results_dir}{variant}.pascal.gz --outdir={results_dir} --customdir=./resources/1kg_AFR/ --custom=AFR --runpathway=on --genescoring=max > {results_dir}{variant}.log"))
  }else if(gene_scoring == 'sum'){
    system(glue::glue("./Pascal --pval {results_dir}{variant}.pascal.gz --outdir={results_dir} --customdir=./resources/1kg_AFR/ --custom=AFR --runpathway=on --genescoring=sum  > {results_dir}{variant}.log"))
  }else{
    stop('gene_scoring misspecified')
  }
  setwd(cur_wd)

}

args <- commandArgs(trailingOnly = T)
gwas_file <- args[[1]]
variant <- args[[2]]
bim_1KG <- args[[3]]
scratch_dir <- args[[4]]
Map1KG(bim_1KG,'ALL',scratch_dir)
results_dir <- args[[5]]
pascal_path <- args[[6]]
RunPASCAL(gwas_file,variant,scratch_dir,results_dir,pascal_path,filt = T)

