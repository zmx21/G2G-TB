library(glue)
library(dplyr)
#Convert VCF file to AA list for each sample with non-synonymous variants
WriteAATable <- function(cur_file,locus_to_excl,snpEff,cur_id,out_path){
  DATA_DIR <- paste0(strsplit(cur_file,split = '/')[[1]][1:(length(strsplit(cur_file,split = '/')[[1]]) - 1)],collapse = '/')
  cur_prefix <- strsplit(cur_file,split = '/')[[1]][length(strsplit(cur_file,split = '/')[[1]])]
  cur_file_unzip <- gsub(pattern = '.gz',replacement = '',cur_prefix)
  system(glue::glue("mkdir -p {DATA_DIR}/rep_filt/"))
  
  out_dir <- gsub(x=out_path,pattern = paste0(cur_id,'.txt'),replacement = '')
  system(glue::glue("mkdir -p {out_dir}"))
  
  #Loop through all possible reading frames for each nucleotide variant
  cur_annot_table <- list()
  loop_tbl = T
  cnt <- 0
  #Loop through all frames
  while(loop_tbl){
    #Remove SNPs in repetitive regions 
    system(glue::glue("bedtools subtract -header -A -a {cur_file} -b {locus_to_excl} > {DATA_DIR}/rep_filt/{cur_file_unzip}"))
    #Get SNP consequences 
    cur_table <- system(glue::glue("java -jar {snpEff} extractFields {DATA_DIR}/rep_filt/{cur_file_unzip} POS REF ALT ANN[{cnt}].EFFECT ANN[{cnt}].GENE EFF[{cnt}].AA GEN[0].GT GEN[0].FREQ> {DATA_DIR}{cur_file_unzip}.tbl"),intern = T)
    cur_table_parsed <- data.table::fread(glue::glue("{DATA_DIR}{cur_file_unzip}.tbl"),sep = '\t',fill = F)
    colnames(cur_table_parsed) <- c('POS','REF','ALT','EFFECT','GENE','AA_Change','GENOTYPE','FREQ')
    cur_table_parsed$FREQ <- sapply(cur_table_parsed$FREQ,function(x) as.numeric(gsub(x=x,pattern = '%',replacement = '')))
    #End loop if there are no longer have any frames to search on. 
    if(!all(is.na(cur_table_parsed$EFFECT))){
      cnt <- cnt + 1
      cur_annot_table <- c(cur_annot_table,list(cur_table_parsed))
    }else{
      loop_tbl <- F
    }
  }
  #Cleanup
  system(glue::glue("rm {DATA_DIR}{cur_file_unzip}.tbl"))
  if(length(cur_annot_table) == 0){
    system(glue::glue("touch {out_path}"))
    return(NA)
  }
  cur_annot_table_full <- do.call(rbind,cur_annot_table) %>% dplyr::filter(EFFECT != 'synonymous_variant' & AA_Change != '')
  data.table::fwrite(cur_annot_table_full,file = out_path)
}

#Convert VCF file to Nucleotide list for each sample with synonymous variants
WriteNucTableSyn <- function(cur_file,locus_to_excl,snpEff,cur_id,out_path){
  DATA_DIR <- paste0(strsplit(cur_file,split = '/')[[1]][1:(length(strsplit(cur_file,split = '/')[[1]]) - 1)],collapse = '/')
  cur_prefix <- strsplit(cur_file,split = '/')[[1]][length(strsplit(cur_file,split = '/')[[1]])]
  cur_file_unzip <- gsub(pattern = '.gz',replacement = '',cur_prefix)
  
  system(glue::glue("mkdir -p {DATA_DIR}/rep_filt/"))
  
  out_dir <- gsub(x=out_path,pattern = paste0(cur_id,'.txt'),replacement = '')
  system(glue::glue("mkdir -p {out_dir}"))
  
  #Loop through all possible reading frames for each nucleotide variant
  cur_annot_table <- list()
  loop_tbl = T
  cnt <- 0
  #Loop through all frames
  while(loop_tbl){
    #Remove SNPs in repetitive regions 
    system(glue::glue("bedtools subtract -header -A -a {cur_file} -b {locus_to_excl} > {DATA_DIR}/rep_filt/{cur_file_unzip}"))
    #Get SNP consequences 
    cur_table <- system(glue::glue("java -jar {snpEff} extractFields {DATA_DIR}/rep_filt/{cur_file_unzip} POS REF ALT ANN[{cnt}].EFFECT ANN[{cnt}].GENE GEN[0].GT GEN[0].FREQ > {DATA_DIR}{cur_file_unzip}.tbl"),intern = T)
    cur_table_parsed <- data.table::fread(glue::glue("{DATA_DIR}{cur_file_unzip}.tbl"),sep = '\t',fill = F)
    colnames(cur_table_parsed) <- c('POS','REF','ALT','EFFECT','GENE','GENOTYPE','FREQ')
    cur_table_parsed$FREQ <- sapply(cur_table_parsed$FREQ,function(x) as.numeric(gsub(x=x,pattern = '%',replacement = '')))
    
    #End loop if there are no longer have any frames to search on. 
    if(!all(is.na(cur_table_parsed$EFFECT))){
      cnt <- cnt + 1
      cur_annot_table <- c(cur_annot_table,list(cur_table_parsed))
    }else{
      loop_tbl <- F
    }
  }
  #Cleanup
  system(glue::glue("rm {DATA_DIR}{cur_file_unzip}.tbl"))
  if(length(cur_annot_table) == 0){
    system(glue::glue("touch {out_path}"))
    return(NA)
  }
  cur_annot_table_full <- do.call(rbind,cur_annot_table)
  syn_variants <- dplyr::filter(cur_annot_table_full,EFFECT == 'synonymous_variant')
  data.table::fwrite(syn_variants,file = out_path)
}


args <- commandArgs(trailingOnly = TRUE) 
vcf_file <- args[[1]]
locus_to_excl <- args[[2]]
snpEff <- args[[3]]
cur_id <- args[[4]]
out_path <- args[[5]]
out_path_syn <- args[[6]]

#Get VCF Files and write into AA tables
WriteAATable(vcf_file,locus_to_excl,snpEff,cur_id,out_path)
WriteNucTableSyn(vcf_file,locus_to_excl,snpEff,cur_id,out_path_syn)
