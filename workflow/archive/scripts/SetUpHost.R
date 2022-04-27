library(dplyr)
SetUpHost <- function(Genotyping_DIR,Pheno_DIR,VCF_Path,Out_Path,excl_regions,imputed = T,SetUpPheno = T,maf_thresh=0.05,covars_discrete_to_incl = c('Patient_Sex','HIV_Status'),covars_numeric_to_incl = c('Age'),n_PC = 3,not_chr_6 = F,col_names=F,X_Chr = T,n_cores){
  system(glue::glue("mkdir -p {Out_Path}"))
  ### Parameters ###
  hwe_thresh = 1e-6 #Hardy-Weinberg Equilibirium
  king_degree <- '2' #Kingship degree (relatives to exclude)
  missing <- 0.1
  KING <- '../software/king'
  
  if(X_Chr){
    X_VCF_Path <- gsub(x=VCF_Path,pattern = '.vcf.gz',replacement = '.chrX.vcf.gz')
  }
  
  
  # ### SNP Filtering ###
  # if(imputed & !file.exists(glue::glue("{VCF_Path}.info"))){
  #   #Write out info score
  #   system(
  #     glue::glue(
  #       "{bcftools} query -f '%INFO/INFO\n' {VCF_Path} > {VCF_Path}.info"
  #     )
  #   )
  #   #Write out MAF
  #   system(
  #     glue::glue(
  #       "{bcftools} query -f '%INFO/RefPanelAF\n' {VCF_Path} > {VCF_Path}.AF"
  #     )
  #   )
  #   if(X_Chr){
  #     system(
  #       glue::glue(
  #         "{bcftools} query -f '%INFO/INFO\n' {X_VCF_Path} > {X_VCF_Path}.info"
  #       )
  #     )
  #     system(
  #       glue::glue(
  #         "{bcftools} query -f '%INFO/RefPanelAF\n' {X_VCF_Path} > {X_VCF_Path}.AF"
  #       )
  #     )
  #   }
  # }
  # 
  # recode_vcf_path <- gsub(x=VCF_Path,pattern = '.vcf.gz',replacement = '.recode')
  # 
  # if(!file.exists(paste0(recode_vcf_path,'.bim'))){
  #   #Set missing ID
  #   system(
  #     glue::glue(
  #       "plink2 --vcf {VCF_Path} --threads {n_cores} --keep-allele-order --make-bed --set-missing-var-ids @:#[b37]\\$r,\\$a --const-fid --out {recode_vcf_path}"
  #     )
  #   )
  #   
  #   if(X_Chr){
  #     recode_vcf_path_X <- gsub(x=X_VCF_Path,pattern = '.vcf.gz',replacement = '.recode')
  #     
  #     system(
  #       glue::glue(
  #         "plink2 --vcf {X_VCF_Path} --threads {n_cores} --keep-allele-order --make-bed --set-missing-var-ids @:#[b37]\\$r,\\$a --const-fid --out {recode_vcf_path_X}.tmp"
  #       )
  #     )
  #     system(
  #       glue::glue(
  #         "plink --bfile {recode_vcf_path_X}.tmp --threads {n_cores} --keep-allele-order --make-bed --impute-sex --out {recode_vcf_path_X}.tmp2"
  #       )
  #     )
  #     # system(glue::glue("rm {recode_vcf_path_X}.tmp.*"))
  #   }
  # }
  # #Calculate MAF
  # system(
  #   glue::glue(
  #     "plink2 --bfile {recode_vcf_path} --freq --threads {n_cores} --out {Out_Path}TB_DAR_Imputed"
  #   )
  # )
  # 
  # 
  # #For X Chr, Split into Males and Females
  # if(X_Chr){
  #   recode_vcf_path_X <- gsub(x=X_VCF_Path,pattern = '.vcf.gz',replacement = '.recode')
  #   
  #   X_fam_file <- data.table::fread(glue::glue('{recode_vcf_path_X}.tmp2.fam'))    
  #   males <-  dplyr::filter(X_fam_file,V5 == 1)
  #   females <-  dplyr::filter(X_fam_file,V5 == 2)
  #   data.table::fwrite(males %>% dplyr::select(V1,V2),file = glue::glue('{recode_vcf_path_X}.males.txt'),col.names = F,row.names = F,quote = F,sep = ' ')
  #   data.table::fwrite(females %>% dplyr::select(V1,V2),file = glue::glue('{recode_vcf_path_X}.females.txt'),col.names = F,row.names = F,quote = F,sep = ' ')
  #   
  #   system(
  #     glue::glue(
  #       "plink2 --bfile {recode_vcf_path_X}.tmp2 --threads {n_cores} --keep {recode_vcf_path_X}.males.txt --make-bed --out {recode_vcf_path_X}.males.tmp"
  #     )
  #   )
  #   
  #   system(
  #     glue::glue(
  #       "plink2 --bfile {recode_vcf_path_X}.tmp2 --threads {n_cores} --keep {recode_vcf_path_X}.females.txt --make-bed --out {recode_vcf_path_X}.females.tmp"
  #     )
  #   )
  #   #MAF,missingness,and HWE filter for each sex (no HWE in males)
  #   system(
  #     glue::glue(
  #       "plink2 --bfile {recode_vcf_path_X}.males.tmp --threads {n_cores} --maf {maf_thresh} --geno {missing} --make-bed --keep-allele-order --out {recode_vcf_path_X}.males.tmp2"
  #     )
  #   )
  #   
  #   system(
  #     glue::glue(
  #       "plink2 --bfile {recode_vcf_path_X}.females.tmp --threads {n_cores} --hwe {hwe_thresh} --maf {maf_thresh} --geno {missing} --make-bed --keep-allele-order --out {recode_vcf_path_X}.females.tmp2"
  #     )
  #   )
  #   
  #   male_variants <- data.table::fread(glue::glue("{recode_vcf_path_X}.males.tmp2.bim"))    
  #   female_variants <- data.table::fread(glue::glue("{recode_vcf_path_X}.females.tmp2.bim"))    
  #   consensus_variants <- intersect(male_variants$V2,female_variants$V2)
  #   write(consensus_variants,glue::glue("{recode_vcf_path_X}.consensus"))
  #   
  #   system(
  #     glue::glue(
  #       "plink2 --bfile {recode_vcf_path_X}.males.tmp2 --threads {n_cores} --extract {recode_vcf_path_X}.consensus --make-bed --keep-allele-order --out {recode_vcf_path_X}.males"
  #     )
  #   )
  #   system(
  #     glue::glue(
  #       "plink2 --bfile {recode_vcf_path_X}.females.tmp2 --threads {n_cores} --extract {recode_vcf_path_X}.consensus --make-bed --keep-allele-order --out {recode_vcf_path_X}.females"
  #     )
  #   )
  #   
  #   system(
  #     glue::glue(
  #       "plink --bfile {recode_vcf_path_X}.males --bmerge {recode_vcf_path_X}.females --threads {n_cores} --indiv-sort f {recode_vcf_path}.fam --make-bed --keep-allele-order --out {recode_vcf_path_X}"
  #     )
  #   )
  #   # system(glue::glue("rm {recode_vcf_path_X}.*tmp*"))
  #   
  # }
  # 
  # 
  # #Write out snps to excl
  # if(imputed){
  #   tmp_bim_file <- data.table::fread(glue::glue("{recode_vcf_path}.bim"))
  #   info_file_header <- data.table::fread(glue::glue("{VCF_Path}.info"),nrows = 1)
  #   if('Rsq' %in% colnames(info_file_header)){
  #     info_file <- data.table::fread(glue::glue("{VCF_Path}.info"),select = 'Rsq')
  #     snps_to_excl <- tmp_bim_file$V2[which(info_file$Rsq < 0.8)]
  #   }else{
  #     info_file <- data.table::fread(glue::glue("{VCF_Path}.info"))
  #     snps_to_excl <- tmp_bim_file$V2[which(info_file$V1 < 0.8)]
  #   }
  #   write(snps_to_excl,file = glue::glue("{Out_Path}to_excl.txt"))
  #   system(
  #     glue::glue(
  #       "plink2 --bfile {recode_vcf_path} --threads {n_cores} --exclude {Out_Path}to_excl.txt --make-bed --keep-allele-order --out {Out_Path}TB_DAR_Imputed_info_filt"
  #     )
  #   )
  #   #Cleanup
  #   system(glue::glue('rm {Out_Path}to_excl.txt'))
  #   
  #   if(X_Chr){
  #     tmp_bim_file <- data.table::fread(glue::glue("{recode_vcf_path_X}.bim"))
  #     info_file_header <- data.table::fread(glue::glue("{X_VCF_Path}.info"),nrows = 1)
  #     if('Rsq' %in% colnames(info_file_header)){
  #       info_file <- data.table::fread(glue::glue("{X_VCF_Path}.info"),select = 'Rsq')
  #       snps_to_excl <- tmp_bim_file$V2[which(info_file$Rsq < 0.8)]
  #     }else{
  #       info_file <- data.table::fread(glue::glue("{X_VCF_Path}.info"))
  #       snps_to_excl <- tmp_bim_file$V2[which(info_file$V1 < 0.8)]
  #     }
  #     write(snps_to_excl,file = glue::glue("{Out_Path}to_excl_X.txt"))
  #     system(
  #       glue::glue(
  #         "plink2 --bfile {recode_vcf_path_X} --threads {n_cores} --exclude {Out_Path}to_excl_X.txt --make-bed --keep-allele-order --out {Out_Path}TB_DAR_Imputed_chrX_info_maf_filt"
  #       )
  #     )
  #     #Cleanup
  #     system(glue::glue('rm {Out_Path}to_excl_X.txt'))
  #     # system(glue::glue('rm {Out_Path}*males*'))
  #     
  #   }
  #   
  #   
  # }else{
  #   system(
  #     glue::glue(
  #       "plink2 --bfile {recode_vcf_path} --threads {n_cores} --make-bed --keep-allele-order --out {Out_Path}TB_DAR_Imputed_info_filt"
  #     )
  #   )
  # }
  # 
  # #MAF and HWE filter on Imputed file
  # system(
  #   glue::glue(
  #     "plink2 --bfile {Out_Path}TB_DAR_Imputed_info_filt --threads {n_cores} --hwe {hwe_thresh} --maf {maf_thresh} --make-bed --keep-allele-order --out {Out_Path}TB_DAR_Imputed_info_maf_filt"
  #   )
  # )
  # 
  # #### Kingship ####
  # system(
  #   glue::glue(
  #     "plink --bfile {Out_Path}TB_DAR_Imputed_info_maf_filt --threads {n_cores} --make-bed --out {Out_Path}tmp"
  #   )
  # )
  # fam_file <- data.table::fread(glue::glue("{Out_Path}tmp.fam"))
  # fam_file$V1 <- fam_file$V2
  # data.table::fwrite(fam_file,sep = ' ',quote = F,col.names = F,row.names = F,file = glue::glue("{Out_Path}tmp.fam"))
  # 
  # system(glue::glue("{KING} --related --degree {king_degree} -b {Out_Path}tmp.bed --prefix {Out_Path}king.related"))
  
  ### Set-up Phenotype ###
  if(SetUpPheno){
    #Load phenotype file
    metadata <- readxl::read_xlsx(glue::glue("{Pheno_DIR}metadata_2020-10-16_17-36-52.xlsx"))
    x_ray_scores <- data.table::fread(glue::glue("{Pheno_DIR}metadata_combined_genomes_022021_inclXray.txt"))
    x_ray_scores$PATIENT_ID <- as.character(x_ray_scores$PATIENT_ID)
    x_ray_scores$Xray_score <- GWAS.utils::trans_inv_normal(x_ray_scores$Xray_score)
    
    colnames(metadata) <- sapply(colnames(metadata),function(x) gsub(pattern = ' ',replacement = '_',x=x))
    
    #Extract relevant columns
    pheno <- dplyr::select(metadata,-G_NUMBER,-Patient_Birthdate) %>% dplyr::left_join(x_ray_scores %>% dplyr::select(PATIENT_ID,Xray_score))
    
    #Transform variables 
    pheno$Age <- as.numeric(pheno$Age)
    pheno$BMI <- as.numeric(pheno$BMI)
    pheno$BMI <- cut(pheno$BMI, c(min(pheno$BMI),18,25,max(pheno$BMI)),labels = c('Underweight','Normal','Overweight'))
    pheno$TB_Score <- as.numeric(pheno$TB_Score)
    pheno$Patient_Household_Size <- as.numeric(pheno$Patient_Household_Size)
    pheno[,grepl(pattern = 'Symptoms_Duration',colnames(pheno))] <- apply(pheno[,grepl(pattern = 'Symptoms_Duration',colnames(pheno))],2,as.numeric)
    pheno[,grepl(pattern = 'Measurements',colnames(pheno))] <- apply(pheno[,grepl(pattern = 'Measurements',colnames(pheno))],2,as.numeric)
    
    #Encode unknown as NA
    pheno$TB_RF_Smoking[pheno$TB_RF_Smoking == 'unknown'] <- NA
    pheno$TB_More1_TB_Site[pheno$TB_More1_TB_Site == 'unknown'] <- NA
    
    #binarize categorical variables
    pheno$TB_RF_Smoking <- factor(pheno$TB_RF_Smoking,labels = c(0,1)) 
    pheno$Patient_Sex <- factor(pheno$Patient_Sex,labels = c(0,1)) #female is 0, male is 1
    pheno$Glucose_Fasting_Gluc <- factor(pheno$Glucose_Fasting_Gluc,labels = c(0,1))
    pheno$General_Examination_Malnutrition <- factor(pheno$General_Examination_Malnutrition,labels = c(0,1))
    pheno$Symptoms_Abdominal_Pain <- factor(pheno$Symptoms_Abdominal_Pain,labels = c(0,1))
    pheno$Symptoms_Chest_Pain <- factor(pheno$Symptoms_Chest_Pain,labels = c(0,1))
    pheno$Symptoms_Cough <- factor(pheno$Symptoms_Cough,labels = c(0,1))
    pheno$Symptoms_Fever <- factor(pheno$Symptoms_Fever,labels = c(0,1))
    pheno$Symptoms_Hemoptysis <- factor(pheno$Symptoms_Hemoptysis,labels = c(0,1))
    pheno$Symptoms_Night_Sweat <- factor(pheno$Symptoms_Night_Sweat,labels = c(0,1))
    pheno$Symptoms_Weight_Loss <- factor(pheno$Symptoms_Weight_Loss,labels = c(0,1))
    
    
    #Check which genotyped samples have pheno information
    fam_file_genotyped <- data.table::fread(glue::glue('{Genotyping_DIR}TB_DAR_Genotyping/PLINK_QCed_Autosomes_Top/BED/Fellay_0620.fam'))
    genotyped_IDs <- sapply(fam_file_genotyped$V2,function(x) ifelse(grepl(x=x,pattern = '-'),strsplit(x=x,split = '-')[[1]][1],x))
    
    sequenced_ids <- data.table::fread(glue::glue("{Genotyping_DIR}training_split.txt"),header = F)$V1
    sequenced_ids <- c(sequenced_ids,data.table::fread(glue::glue("{Genotyping_DIR}testing_split.txt"),header = F)$V1)
    sequenced_ids <- sapply(sequenced_ids,function(x) strsplit(x=x,split = '\\.')[[1]][2])
    
    print(paste0('Samples with Pheno:',length(intersect(pheno$PATIENT_ID,c(genotyped_IDs,sequenced_ids)))))
  }
  ### Sample Filtering ###
  #Exclude duplicate samples and 1 indiv. from the related sample pairs
  king_related <- data.table::fread(glue::glue("{Out_Path}king.related.kin0")) %>% dplyr::select(ID1,ID2,Kinship,InfType)
  king_related$Simple_ID1 <- sapply(king_related$ID1,function(x){ 
    if(grepl(x=x,pattern = '-')){
      return(strsplit(x=strsplit(x=x,split = '-')[[1]][1],split = '_')[[1]][2])
    }else if(grepl(x=x,pattern = 'WGS_Fellay.')){
      return(strsplit(x=x,split = 'WGS_Fellay.')[[1]][2])
    }else{
      return(strsplit(x=x,split = '_')[[1]][2])
    }
  })
  king_related$Simple_ID2 <- sapply(king_related$ID2,function(x){ 
    if(grepl(x=x,pattern = '-')){
      return(strsplit(x=strsplit(x=x,split = '-')[[1]][1],split = '_')[[1]][2])
    }else if(grepl(x=x,pattern = 'WGS_Fellay.')){
      return(strsplit(x=x,split = 'WGS_Fellay.')[[1]][2])
    }else{
      return(strsplit(x=x,split = '_')[[1]][2])
    }
  })
  pheno_present_ID1 <- king_related$Simple_ID1 %in% pheno$PATIENT_ID
  pheno_present_ID2 <- king_related$Simple_ID2 %in% pheno$PATIENT_ID
  king_related$Pheno_Avail <- mapply(function(x,y) paste(x,y),pheno_present_ID1,pheno_present_ID2)
  data.table::fwrite(king_related %>% dplyr::select(ID1=Simple_ID1,ID2=Simple_ID2,Kinship,InfType,Pheno_Avail),file = glue::glue("{Out_Path}relatedness.csv"),sep = ',')
  
  #Remove Dup samples
  dup_samples <- dplyr::filter(king_related,InfType == 'Dup/MZ')
  king_excl <- c(dup_samples$ID1,dup_samples$ID2)
  
  #keep 1 of the samples for relatives (the one with phenotype if the other one doesn't have pheno)
  king_related$Sample_To_Exclu <- mapply(FUN = function(x,y){
    if(sum(c(x,y)) == 2){
      return(1)
    }else if (sum(c(x,y)) == 0){
      return(c(1,2))
    }else{
      return(which(c(x,y)))
    }
  },pheno_present_ID1,pheno_present_ID2)
  relatives <- dplyr::filter(king_related,InfType != 'Dup/MZ')
  relatives_to_excl <- unlist(lapply(1:nrow(relatives),function(i) c(relatives$ID1[i],relatives$ID2[i])[relatives$Sample_To_Exclu[[i]]]))
  king_excl <- c(king_excl,relatives_to_excl)
  
  #PCA outliers to remove (WARNING: Checked manually beforehand on which PC1 cutoff to use)
  source('./scripts/Run_PCA.R')
  #Prepare 1KG sample table
  PCA_Out_Path <- paste0(Out_Path,'PCA')

  RunPCA1KG(VCF_Path,PCA_Out_Path)
  
  super_pop <- data.table::fread('~/G2G_TB/data/1000_Genomes/20131219.populations.tsv')
  tbl_1KG <- data.table::fread('~/G2G_TB/data/1000_Genomes/20130606_g1k.ped') %>% dplyr::left_join(super_pop,by = c('Population' = 'Population Code'))
  
  #Load PCA results
  pc_df <- data.table::fread(glue::glue("{PCA_Out_Path}_1KG_TBDAR_consensus_pruned.eigenvec"))
  colnames(pc_df) <- c('FID','IID',paste0('PC',seq(1,ncol(pc_df) - 2)))
  #Remove ID Prefix
  pc_df$IID <- sapply(pc_df$IID,function(x) ifelse(grepl(pattern='_',x=x),strsplit(x=x,split = '_')[[1]][2],x))
  #Plot with all Superpopulations
  merged_df <- pc_df %>%
    dplyr::left_join(tbl_1KG,c('IID'='Individual ID'))
  merged_df$`Data Set` <- as.factor(ifelse(is.na(merged_df$Population),'TB DAR','1000 Genomes'))
  merged_df$`Super Population` <- as.factor(as.character(merged_df$`Super Population`))
  merged_df$`Super Population` <-factor(merged_df$`Super Population` ,
                                        levels = rev(levels(merged_df$`Super Population`)))
  
  eigen_val <- as.vector(data.table::fread(glue::glue("{PCA_Out_Path}_1KG_TBDAR_consensus_pruned.eigenval")))
  var_explained <- eigen_val$V1 / sum(eigen_val$V1) * 100
  merged_df$simple_ID <- sapply(merged_df$IID,function(x) ifelse(grepl(x=x,pattern = '-'),strsplit(x=x,split = '-')[[1]][1],x))
  pc_plot_thousand_genome_PC1_PC2 <- ggplot(data = merged_df) +
    aes(x=PC1,y=PC2,color=`Super Population`,shape=`Data Set`) +
    geom_point() + scale_shape_manual(values = c(20,4)) +
    xlab(paste0('PC1 (',signif(var_explained[1],2),'%)'))  +
    ylab(paste0('PC2 (',signif(var_explained[2],2),'%)')) +
    geom_text_repel(data= dplyr::filter(merged_df,`Data Set` == 'TB DAR' & PC1 < 0),aes(label=simple_ID))
  pc_plot_thousand_genome_PC3_PC4 <- ggplot(data = merged_df) +
    aes(x=PC3,y=PC4,color=`Super Population`,shape=`Data Set`) +
    geom_point() + scale_shape_manual(values = c(20,4)) +
    xlab(paste0('PC3 (',signif(var_explained[3],2),'%)'))  +
    ylab(paste0('PC4 (',signif(var_explained[4],2),'%)')) #+
  pc_plot_thousand_genome <- ggpubr::ggarrange(plotlist = list(pc_plot_thousand_genome_PC1_PC2,pc_plot_thousand_genome_PC3_PC4),ncol = 2)
  
  tb_dar_PCA <- dplyr::filter(merged_df,`Data Set` == 'TB DAR')
  pca_to_excl <- dplyr::filter(tb_dar_PCA,PC1 < 0)$IID
  
  if(X_Chr){
    #Check Sex Mismatch 
    imputed_fam_file_X <- data.table::fread(glue::glue("{Out_Path}TB_DAR_Imputed_chrX_info_maf_filt.fam"))
    imputed_IDs <- sapply(imputed_fam_file_X$V2,function(x) strsplit(x=x,split = '_')[[1]][2])
    imputed_IDs <- sapply(imputed_IDs,function(x) ifelse(grepl(x=x,pattern = '-'),strsplit(x=x,split = '-')[[1]][1],x))
    
    metadata_sex <- metadata$Patient_Sex[match(imputed_IDs,metadata$PATIENT_ID)]
    imputed_sex <- sapply(imputed_fam_file_X$V5,function(x){
      if(x==1){
        'male'
      }else if(x==2){
        'female'
      }else{
        NA
      }
    })
    imputed_fam_file <- data.table::fread(glue::glue("{Out_Path}TB_DAR_Imputed_info_maf_filt.fam"))
    sex_mismatch <- imputed_fam_file_X$V2[which(metadata_sex != imputed_sex)]
    undefined_sex <- setdiff(imputed_fam_file$V2,imputed_fam_file_X$V2)
    fam_file_to_excl <- sapply(imputed_fam_file$V2,function(x) x %in% king_excl | x%in% sex_mismatch | x %in% undefined_sex | any(sapply(pca_to_excl,function(y) grepl(x=x,pattern = y))))
    
  }else{
    imputed_fam_file <- data.table::fread(glue::glue("{Out_Path}TB_DAR_Imputed_info_maf_filt.fam"))
    fam_file_to_excl <- sapply(imputed_fam_file$V2,function(x) x %in% king_excl | any(sapply(pca_to_excl,function(y) grepl(x=x,pattern = y))))
    
  }
  
  #Write out filtered genotype file for GWAS
  data.table::fwrite(imputed_fam_file[fam_file_to_excl,c(1,2)],sep = ' ',col.names = F,row.names = F,file = glue::glue("{Out_Path}samples_to_excl.txt"))
  system(
    glue::glue(
      "plink2 --bfile {Out_Path}TB_DAR_Imputed_info_maf_filt --threads {n_cores} --remove {Out_Path}samples_to_excl.txt --make-bed --keep-allele-order --rm-dup force-first --out {Out_Path}TB_DAR_Imputed_GWAS"
    )
  )
  if(X_Chr){
    system(
      glue::glue(
        "plink2 --bfile {Out_Path}TB_DAR_Imputed_chrX_info_maf_filt --threads {n_cores} --remove {Out_Path}samples_to_excl.txt --make-bed --keep-allele-order --rm-dup force-first --out {Out_Path}TB_DAR_Imputed_chrX_GWAS"
      )
    )
  }
  
  #File for PCA (Prune and remove long-range LD)
  system(
    glue::glue(
      "plink2 --bfile {Out_Path}TB_DAR_Imputed_GWAS --threads {n_cores} --keep-allele-order --new-id-max-allele-len 50 --set-all-var-ids @:#[b37]\\$r,\\$a --rm-dup force-first --make-bed --out {Out_Path}TB_DAR_Imputed_GWAS_PCA_tmp"
    )
  )
  
  system(
    glue::glue(
      "plink2 --bfile {Out_Path}TB_DAR_Imputed_GWAS_PCA_tmp --threads {n_cores} --keep-allele-order --indep-pairwise 200 100 0.2 --out {Out_Path}TB_DAR_Imputed_GWAS_PCA_tmp"
    )
  )
  
  if(not_chr_6){
    system(
      glue::glue(
        "plink2 --bfile {Out_Path}TB_DAR_Imputed_GWAS_PCA_tmp --threads {n_cores} --exclude range {excl_regions} --not-chr 6 --keep-allele-order --extract {Out_Path}TB_DAR_Imputed_GWAS_PCA_tmp.prune.in --make-bed --out {Out_Path}TB_DAR_Imputed_GWAS_PCA"
      )
    )
    
  }else{
    system(
      glue::glue(
        "plink2 --bfile {Out_Path}TB_DAR_Imputed_GWAS_PCA_tmp --threads {n_cores} --exclude range {excl_regions} --keep-allele-order --extract {Out_Path}TB_DAR_Imputed_GWAS_PCA_tmp.prune.in --make-bed --out {Out_Path}TB_DAR_Imputed_GWAS_PCA"
      )
    )
    
  }
  
  #Calculate PCs 
  system(
    glue::glue(
      "gcta64 --bfile {Out_Path}TB_DAR_Imputed_GWAS --make-grm --out {Out_Path}TB_DAR_Imputed_GWAS --thread-num {n_cores}"
    )
  )
  
  system(
    glue::glue(
      "gcta64 --bfile {Out_Path}TB_DAR_Imputed_GWAS_PCA --make-grm --out {Out_Path}TB_DAR_Imputed_GWAS_PCA --thread-num {n_cores}"
    )
  )
  system(
    glue::glue(
      "gcta64 --grm {Out_Path}TB_DAR_Imputed_GWAS_PCA --pca 20 --out {Out_Path}TB_DAR_Imputed_GWAS_PCA --thread-num {n_cores}"
    )
  )
  
  if(SetUpPheno){
    #Parse sample ID
    samples_for_GWAS <- as.vector(t(imputed_fam_file[!fam_file_to_excl,2]))
    samples_for_GWAS_simple <- sapply(samples_for_GWAS,function(x) strsplit(x=x,split = '_')[[1]][2])
    samples_for_GWAS_simple <- sapply(samples_for_GWAS_simple,function(x) ifelse(grepl(x=x,pattern = '-'),strsplit(x=x,split = '-')[[1]][1],x))
    
    #Join with pheno table, extract relevant covars
    tb_score <- dplyr::left_join(data.frame(ID = samples_for_GWAS_simple,stringsAsFactors = F),pheno %>% dplyr::select(PATIENT_ID,TB_Score),by=c('ID'='PATIENT_ID'))$TB_Score
    xray_score <- dplyr::left_join(data.frame(ID = samples_for_GWAS_simple,stringsAsFactors = F),pheno %>% dplyr::select(PATIENT_ID,Xray_score),by=c('ID'='PATIENT_ID'))$Xray_score
    
    covars_discrete <- dplyr::left_join(data.frame(ID = samples_for_GWAS_simple,stringsAsFactors = F),pheno %>% 
                                          dplyr::select(c('PATIENT_ID',covars_discrete_to_incl)),by=c('ID'='PATIENT_ID')) %>% dplyr::select(-ID)
    covars_numeric <- dplyr::left_join(data.frame(ID = samples_for_GWAS_simple,stringsAsFactors = F),pheno %>% 
                                         dplyr::select(c('PATIENT_ID',covars_numeric_to_incl)),by=c('ID'='PATIENT_ID')) %>% dplyr::select(-ID)
    PCs <- data.table::fread(glue::glue("{Out_Path}TB_DAR_Imputed_GWAS_PCA.eigenvec"))
    colnames(PCs) <- c('FID','IID',paste0('PC',seq(1,ncol(PCs) - 2)))
    
    #Top 3 PCs and age
    data.table::fwrite(cbind(PCs[,1:(2+n_PC)],covars_numeric),sep = ' ',col.names = col_names,na = 'NA',file = glue::glue("{Out_Path}covars_numeric"),quote = F)
    
    #sex, BMI and HIV status
    data.table::fwrite(cbind(PCs[,1:2],covars_discrete),sep = ' ',col.names = col_names,na = 'NA',file = glue::glue("{Out_Path}covars_discrete"),quote = F)
    
    #pheno 
    data.table::fwrite(cbind(PCs[,1:2],data.frame(tb_score=tb_score)),sep = ' ',col.names = col_names,na = 'NA',file = glue::glue("{Out_Path}tb_score"),quote = F)
    
    #x-ray score 
    data.table::fwrite(cbind(PCs[,1:2],data.frame(xray_score=xray_score)),sep = ' ',col.names = col_names,na = 'NA',file = glue::glue("{Out_Path}xray_score"),quote = F)
    
    
    
    symptoms <- dplyr::left_join(data.frame(ID = samples_for_GWAS_simple,stringsAsFactors = F),pheno %>% dplyr::select(PATIENT_ID,Symptoms_Cough,Symptoms_Fever,Symptoms_Weight_Loss,Symptoms_Abdominal_Pain,Symptoms_Hemoptysis,Symptoms_Chest_Pain,Symptoms_Night_Sweat),by=c('ID'='PATIENT_ID'))
    symptoms_num <- apply(symptoms[,-1],1,function(x) sum(as.numeric(as.character(x))))
    data.table::fwrite(cbind(PCs[,1:2],data.frame(symptoms_num=symptoms_num)),sep = ' ',col.names = col_names,na = 'NA',file = glue::glue("{Out_Path}symptoms_num"),quote = F)
    # summary(lm('symptoms_num ~ PC1 + PC2 + PC3 + hiv_status + age + sex',cbind(data.frame(symptoms_num=symptoms_num,
    #                                                                   sex=covars_discrete$Patient_Sex,age = covars_numeric$Age,
    #                                                                   hiv_status = covars_discrete$HIV_Status),PCs[,-c(1,2)])))
  }
  # system(glue::glue('rm {Out_Path}*tmp*'))
  if(X_Chr){
    system(
      glue::glue(
        "plink --bfile {Out_Path}TB_DAR_Imputed_GWAS --bmerge {Out_Path}TB_DAR_Imputed_chrX_GWAS --make-bed --out {Out_Path}TB_DAR_Imputed_GWAS_Full"
      )
    )
    
    sapply(c('bim','bed','fam'),function(x)system(
      glue::glue(
        "mv {Out_Path}TB_DAR_Imputed_GWAS_Full.{x} {Out_Path}TB_DAR_Imputed_GWAS.{x}"
      )
    ))
    
  }
}

SetUpHLA <- function(VCF_File,Out_Path,MAF_Thresh = 0.01,Info_Thresh = 0.8){
  system(glue::glue("mkdir -p {Out_Path}"))
  info_file <- data.table::fread(gsub(VCF_File,pattern = '.dose.vcf.gz',replacement = '.info'))
  
  hla_alleles <- dplyr::filter(info_file,Rsq > Info_Thresh & grepl('HLA',SNP)) %>% dplyr::select(SNP)
  hla_AA <- dplyr::filter(info_file,Rsq > Info_Thresh & grepl('AA',SNP)) %>% dplyr::select(SNP)
  
  write(hla_AA$SNP,file = glue::glue("{Out_Path}AA_to_keep.txt"))
  write(hla_alleles$SNP,file = glue::glue("{Out_Path}Alleles_to_keep.txt"))
  
  system(glue::glue("plink2 --vcf {VCF_File} --extract {Out_Path}AA_to_keep.txt --maf {MAF_Thresh} --make-bed --out {Out_Path}TB_DAR_HLA_AA"))
  system(glue::glue("plink2 --vcf {VCF_File} --extract {Out_Path}Alleles_to_keep.txt --maf {MAF_Thresh} --make-bed --out {Out_Path}TB_DAR_HLA_Alleles"))
  
}


args <- commandArgs(trailingOnly = TRUE) 
Genotyping_DIR <- args[[1]]
Pheno_DIR <- args[[2]]
VCF_Path <- args[[3]]
HLA_Path <- args[[4]]
OUT_dir <- args[[5]]
excl_regions <- args[[6]]
Host_MAF <- as.numeric(args[[7]])
X_Chr <- as.logical(args[[8]])
n_cores <- as.numeric(args[[9]])

SetUpHost(Genotyping_DIR,Pheno_DIR,VCF_Path,OUT_dir,excl_regions=excl_regions,maf_thresh=Host_MAF,
          covars_discrete_to_incl = c('Patient_Sex','HIV_Status','TB_RF_Smoking','BMI','Patient_Household_Size'),covars_numeric_to_incl = c('Age'),col_names=T,X_Chr = X_Chr,n_cores=n_cores)
SetUpHLA(VCF_File = HLA_Path,Out_Path = OUT_dir)