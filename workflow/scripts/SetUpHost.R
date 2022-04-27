library(dplyr)
library(tidyr)
library(dummies)
library(data.table)

SetUpHost <- function(Metadata_Path,BFILE_Path,Out_Path,excl_regions,covars_discrete_to_incl = c('patient_sex','HIV_status'),covars_numeric_to_incl = c('age','BMI'),pheno_name = c('TB_score','Xray_score'),n_PC = 3,not_chr_6 = T,n_cores){
  system(glue::glue("mkdir -p {Out_Path}"))
  system(glue::glue("mkdir -p {Out_Path}/Host/"))
  
  #File for PCA (Prune and remove long-range LD)
  system(
    glue::glue(
      "plink2 --bfile {BFILE_Path} --threads {n_cores} --keep-allele-order --new-id-max-allele-len 50 --set-all-var-ids @:#[b37]\\$r,\\$a --rm-dup force-first --make-bed --out {Out_Path}/Host/TB_DAR_GWAS_PCA_tmp"
    )
  )
  
  system(
    glue::glue(
      "plink2 --bfile {Out_Path}/Host/TB_DAR_GWAS_PCA_tmp --threads {n_cores} --keep-allele-order --indep-pairwise 200 100 0.2 --out {Out_Path}/Host/TB_DAR_GWAS_PCA_tmp"
    )
  )
  
  if(not_chr_6){
    system(
      glue::glue(
        "plink2 --bfile {Out_Path}/Host/TB_DAR_GWAS_PCA_tmp --threads {n_cores} --exclude range {excl_regions} --not-chr 6 --keep-allele-order --extract {Out_Path}/Host/TB_DAR_GWAS_PCA_tmp.prune.in --make-bed --out {Out_Path}/Host/TB_DAR_GWAS_PCA"
      )
    )
    
  }else{
    system(
      glue::glue(
        "plink2 --bfile {Out_Path}/Host/TB_DAR_GWAS_PCA_tmp --threads {n_cores} --exclude range {excl_regions} --keep-allele-order --extract {Out_Path}/Host/TB_DAR_GWAS_PCA_tmp.prune.in --make-bed --out {Out_Path}/Host/TB_DAR_GWAS_PCA"
      )
    )
    
  }
  
  ### Calculate PCs ###
  system(
    glue::glue(
      "gcta64 --bfile {BFILE_Path} --make-grm --out {Out_Path}/Host/TB_DAR_GWAS --thread-num {n_cores}"
    )
  )
  
  system(
    glue::glue(
      "gcta64 --bfile {Out_Path}/Host/TB_DAR_GWAS_PCA --make-grm --out {Out_Path}/Host/TB_DAR_GWAS_PCA --thread-num {n_cores}"
    )
  )
  system(
    glue::glue(
      "gcta64 --grm {Out_Path}/Host/TB_DAR_GWAS_PCA --pca 20 --out {Out_Path}/Host/TB_DAR_GWAS_PCA --thread-num {n_cores}"
    )
  )
  
  ### Set-up Phenotype ###
  #Load metadata file
  metadata <- data.table::fread(Metadata_Path,na.strings = c("",'unknown'),stringsAsFactors = T) %>%
    dplyr::select(any_of(c('PATIENT_ID',covars_discrete_to_incl,covars_numeric_to_incl))) %>% dplyr::mutate(PATIENT_ID = as.character(PATIENT_ID)) 
  
  #Parse sample ID
  samples_for_GWAS <- data.table::fread(glue::glue("{BFILE_Path}.fam")) %>% dplyr::select(IID=V2)
  samples_for_GWAS_simple <- sapply(samples_for_GWAS$IID,function(x) gsub(x=x,pattern = 'WGS_Fellay.',replacement = ''))
  samples_for_GWAS_simple <- sapply(samples_for_GWAS_simple,function(x) gsub(x=x,pattern = 'Batch1_',replacement = ''))
  samples_for_GWAS_simple <- sapply(samples_for_GWAS_simple,function(x) gsub(x=x,pattern = 'Batch2_',replacement = ''))
  
  if(is.null(covars_discrete_to_incl)){
    covars_discrete <- NULL
  }else{
    covars_discrete <- dplyr::left_join(data.frame(ID = samples_for_GWAS_simple,stringsAsFactors = F),metadata %>% 
                                          dplyr::select(any_of(c('PATIENT_ID',covars_discrete_to_incl))),by=c('ID'='PATIENT_ID')) %>% dplyr::select(-ID)
  }
  
  if(is.null(covars_numeric_to_incl)){
    covars_numeric <- NULL
  }else{
    covars_numeric <- dplyr::left_join(data.frame(ID = samples_for_GWAS_simple,stringsAsFactors = F),metadata %>% 
                                         dplyr::select(any_of(c('PATIENT_ID',covars_numeric_to_incl))),by=c('ID'='PATIENT_ID')) %>% dplyr::select(-ID)
    
  }
  
  PCs <- data.table::fread(glue::glue("{Out_Path}/Host/TB_DAR_GWAS_PCA.eigenvec"))
  colnames(PCs) <- c('FID','IID',paste0('PC',seq(1,ncol(PCs) - 2)))
  
  #Top N PCs and numeric covars
  data.table::fwrite(cbind(PCs[,1:(2+n_PC)],covars_numeric),sep = '\t',col.names = T,na = 'NA',file = glue::glue("{Out_Path}/Host/covars_numeric"))
  
  #discrete covars
  data.table::fwrite(cbind(PCs[,1:2],covars_discrete),sep = '\t',col.names = T,na = 'NA',file = glue::glue("{Out_Path}/Host/covars_discrete"))
  
  for(i in 1:length(pheno_name)){
    metadata_pheno <- data.table::fread(Metadata_Path,na.strings = c("",'unknown'),stringsAsFactors = T) %>%
      dplyr::select(any_of(c('PATIENT_ID',pheno_name[i]))) %>% dplyr::mutate(PATIENT_ID = as.character(PATIENT_ID)) %>% 
      tidyr::drop_na()
    
    #Join with pheno table, extract relevant covars
    pheno_vect <- dplyr::left_join(data.frame(ID = samples_for_GWAS_simple,stringsAsFactors = F),metadata_pheno %>% dplyr::select(PATIENT_ID,any_of(pheno_name[i])),by=c('ID'='PATIENT_ID'))[,pheno_name[i],drop=T]
    
    #pheno 
    if(pheno_name[i] == 'Lineage'){
      pheno_df <- cbind(PCs[,1:2],data.frame(Lineage = pheno_vect)) %>% tidyr::drop_na()
      pheno_df <- cbind(pheno_df[,1:2],dummy.data.frame(pheno_df[,-c(1,2)]))
    }else if(pheno_name[i] == 'Sublineage'){
      pheno_df <- cbind(PCs[,1:2],data.frame(Lineage = pheno_vect)) %>% tidyr::drop_na()
      dummy_data <- dummy.data.frame(pheno_df[,-c(1,2)])
      #Filter Sublineage with less than 30 counts
      cnts <- apply(dummy_data,2,function(x) min(c(sum(x==1),sum(x==0))))
      dummy_data <- dummy_data[,cnts > 30]
      pheno_df <- cbind(pheno_df[,1:2],dummy_data)
    }
    else{
      pheno_df <- cbind(PCs[,1:2],data.frame(pheno=pheno_vect))
      colnames(pheno_df)[3] <- pheno_name[i]
    }
    data.table::fwrite(pheno_df,sep = ' ',col.names = T,na = 'NA',file = glue::glue("{Out_Path}/Host/{pheno_name[i]}"),quote = F)
  }
  
}


SetUpHLA <- function(VCF_File,Out_Path,MAF_Thresh = 0.01,Info_Thresh = 0.8){
  info_file <- data.table::fread(cmd = glue::glue("zcat {gsub(VCF_File,pattern = '.dose.vcf.gz',replacement = '.info.gz')}"))
  
  hla_alleles <- dplyr::filter(info_file,Rsq > Info_Thresh & grepl('HLA',SNP)) %>% dplyr::select(SNP)
  hla_AA <- dplyr::filter(info_file,Rsq > Info_Thresh & grepl('AA',SNP)) %>% dplyr::select(SNP)
  
  write(hla_AA$SNP,file = glue::glue("{Out_Path}Host/AA_to_keep.txt"))
  write(hla_alleles$SNP,file = glue::glue("{Out_Path}Host/Alleles_to_keep.txt"))
  
  system(glue::glue("plink2 --vcf {VCF_File} --extract {Out_Path}Host/AA_to_keep.txt --maf {MAF_Thresh} --make-bed --double-id --out {Out_Path}Host/TB_DAR_HLA_AA"))
  system(glue::glue("plink2 --vcf {VCF_File} --extract {Out_Path}Host/Alleles_to_keep.txt --maf {MAF_Thresh} --make-bed --double-id --out {Out_Path}Host/TB_DAR_HLA_Alleles"))
  
}

args <- commandArgs(trailingOnly = TRUE)
Metadata_Path <- args[[1]]
BFILE_Path <- gsub(args[[2]],pattern = '.bed',replacement = '')
HLA_Path <- args[[3]]
OUT_dir <- args[[4]]
excl_regions <- args[[5]]
n_cores <- as.numeric(args[[6]])

SetUpHost(Metadata_Path,BFILE_Path,OUT_dir,excl_regions=excl_regions,covars_discrete_to_incl = c('patient_sex','HIV_status','TB_RF_smoking'),covars_numeric_to_incl = c('age','BMI'),n_cores=n_cores)
SetUpHLA(VCF_File = HLA_Path,Out_Path = OUT_dir)