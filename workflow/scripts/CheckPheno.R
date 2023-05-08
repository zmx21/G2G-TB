#Credit: Sina Rueger
trans_inv_normal <- function(x, const = 3 / 8) {
  rank.x <- rank(x, na.last = "keep")
  N <- length(na.omit(x))
  qnorm((rank.x - const) / (N - 2 * const + 1))
}


CheckPheno <- function(G2G_Obj,host_snp,AA_variant,lineage,pheno,analysis = 'Both'){
  if(analysis == 'Both'){
    host_genotype <- snpStats::read.plink(bed = paste0('../',G2G_Obj$host_path[lineage]),select.snps = host_snp)
    host_dosage <- as(host_genotype$genotypes, Class = 'numeric')
    host_dosage_ordered <- host_dosage[match(G2G_Obj$both_IDs_to_keep[[lineage]]$FAM_ID,rownames(host_dosage)),]
    
    #Flip allele, REF allele as 0
    host_dosage_ordered <- abs(host_dosage_ordered-2)
    
    #host_dosage_ordered[host_dosage_ordered==2] <- paste0(host_genotype$map$allele.2,host_genotype$map$allele.2)
    #host_dosage_ordered[host_dosage_ordered==1] <- paste0(host_genotype$map$allele.1,host_genotype$map$allele.2)
    #host_dosage_ordered[host_dosage_ordered==0] <- paste0(host_genotype$map$allele.1,host_genotype$map$allele.1)
    
    if(lineage == 'ALL'){
      pathogen_dosage <- G2G_Obj$aa_matrix_full[,AA_variant,drop = F]
      
    }else{
      pathogen_dosage <- G2G_Obj$aa_matrix_filt[[lineage]][,AA_variant,drop = F]
    }
    pathogen_dosage[pathogen_dosage == 2] <- 1
    pathogen_dosage_ordered <- pathogen_dosage[match(G2G_Obj$both_IDs_to_keep[[lineage]]$G_NUMBER,rownames(pathogen_dosage)),,drop = F]
    
    metadata <- data.frame(G_NUMBER = rownames(pathogen_dosage_ordered),Host_SNP = host_dosage_ordered,Pathogen_SNP = as.vector(pathogen_dosage_ordered))
    pathogen_data <- data.table::fread('../../data/pheno/metadata_Sinergia_final_dataset_human_bac_genome_available_QCed.txt')
    pathogen_data$HIV_status[pathogen_data$HIV_status == ''] <- NA
    pathogen_data$TB_RF_smoking[pathogen_data$TB_RF_smoking == 'unknown'] <- NA
    
    pathogen_data_jned <- dplyr::left_join(metadata,pathogen_data)
    Host_PCs <- G2G_Obj$host_PCs$ALL
    Host_PCs$PATIENT_ID <- sapply(Host_PCs$IID,function(x) gsub(x=x,pattern = 'Batch1_',replacement = ''))
    Host_PCs$PATIENT_ID <- sapply(Host_PCs$PATIENT_ID,function(x) gsub(x=x,pattern = 'Batch2_',replacement = ''))
    Host_PCs$PATIENT_ID <- sapply(Host_PCs$PATIENT_ID,function(x) gsub(x=x,pattern = 'WGS_Fellay.',replacement = ''))
    Host_PCs$PATIENT_ID <- as.numeric(Host_PCs$PATIENT_ID)
    
    library(MASS)
    TB_score_data <- na.omit(dplyr::select(pathogen_data_jned,TB_score,age,HIV_status,patient_sex,symptoms_duration_cough_duration,TB_RF_smoking))
    TB_Score_polr <- polr(factor(TB_score) ~., data = TB_score_data,Hess = T)
    TB_Score_step <- stepAIC(TB_Score_polr, direction = "both", trace = FALSE)
    
    Xray_score_data <- na.omit(dplyr::select(pathogen_data_jned,Xray_score,age,HIV_status,patient_sex,symptoms_duration_cough_duration,TB_RF_smoking))
    Xray_score_data$Xray_score_int <- trans_inv_normal(Xray_score_data$Xray_score)
    Xray_score_data <- Xray_score_data %>% dplyr::select(-Xray_score)
    
    Xray_score_lm <- lm(Xray_score_int ~., data = Xray_score_data)
    Xray_score_step <- stepAIC(Xray_score_lm, direction = "both", trace = FALSE)
    
    
    Ct_value_data <- na.omit(dplyr::select(pathogen_data_jned,Ct_value,age,HIV_status,patient_sex,symptoms_duration_cough_duration,TB_RF_smoking))
    Ct_value_data$Ct_value_int <- trans_inv_normal(Ct_value_data$Ct_value)
    Ct_value_data <- Ct_value_data %>% dplyr::select(-Ct_value)
    
    Ct_value_lm <- lm(Ct_value_int ~., data = Ct_value_data)
    Ct_value_step <- stepAIC(Ct_value_lm, direction = "both", trace = FALSE)
    
    
    TB_Score_Data_SNP <- na.omit(dplyr::select(pathogen_data_jned,PATIENT_ID,TB_score,HIV_status,patient_sex,symptoms_duration_cough_duration,TB_RF_smoking,Pathogen_SNP,Host_SNP) %>% dplyr::left_join(Host_PCs))
    
    TB_Score_Data_SNP_Both_None <- dplyr::filter(TB_Score_Data_SNP, (Pathogen_SNP &  Host_SNP) | (!Pathogen_SNP & !Host_SNP))
    TB_Score_Data_SNP_Both_None$Pathogen_Host <- TB_Score_Data_SNP_Both_None$Pathogen_SNP & TB_Score_Data_SNP_Both_None$Host_SNP
    
    
    # print(coeftest(polr('factor(TB_score) ~ HIV_status + symptoms_duration_cough_duration + TB_RF_smoking +  patient_sex + Pathogen_Host', 
    #                     data = TB_Score_Data_SNP_Both_None,Hess = T)))
    
    print(coeftest(polr('factor(TB_score) ~ HIV_status + symptoms_duration_cough_duration + TB_RF_smoking +  patient_sex + Host_SNP * Pathogen_SNP', 
                        data = TB_Score_Data_SNP,Hess = T)))
    print(coeftest(polr('factor(TB_score) ~ HIV_status + symptoms_duration_cough_duration + TB_RF_smoking +  patient_sex + Host_SNP', 
                        data = TB_Score_Data_SNP,Hess = T)))
    print(coeftest(polr('factor(TB_score) ~ HIV_status + symptoms_duration_cough_duration + TB_RF_smoking +  patient_sex + Pathogen_SNP', 
                        data = TB_Score_Data_SNP,Hess = T)))
    
    Xray_score_data_SNP <- na.omit(dplyr::select(pathogen_data_jned,Xray_score,age,HIV_status,patient_sex,symptoms_duration_cough_duration,TB_RF_smoking,Pathogen_SNP,Host_SNP))
    Xray_score_data_SNP$Xray_score_int <- trans_inv_normal(Xray_score_data_SNP$Xray_score)

    Xray_score_data_SNP_Both_None <- dplyr::filter(Xray_score_data_SNP, (Pathogen_SNP &  Host_SNP) | (!Pathogen_SNP & !Host_SNP))
    Xray_score_data_SNP_Both_None$Pathogen_Host <- Xray_score_data_SNP_Both_None$Pathogen_SNP & Xray_score_data_SNP_Both_None$Host_SNP
    
    # print(summary(lm('Xray_score_int ~ age + patient_sex + HIV_status + Pathogen_Host', 
    #                  data = Xray_score_data_SNP_Both_None)))
    
    print(summary(lm('Xray_score_int ~ age + patient_sex + HIV_status + Host_SNP*Pathogen_SNP', 
                     data = Xray_score_data_SNP)))
    print(summary(lm('Xray_score_int ~ age + patient_sex + HIV_status + Host_SNP', 
                     data = Xray_score_data_SNP)))
    print(summary(lm('Xray_score_int ~ age + patient_sex + HIV_status + Pathogen_SNP', 
                     data = Xray_score_data_SNP)))

    Ct_value_data_SNP <- na.omit(dplyr::select(pathogen_data_jned,Ct_value,HIV_status,symptoms_duration_cough_duration,TB_RF_smoking,Pathogen_SNP,Host_SNP))
    Ct_value_data_SNP$Ct_value_int <- trans_inv_normal(Ct_value_data_SNP$Ct_value)
    
    Ct_value_data_SNP_Both_None <- dplyr::filter(Ct_value_data_SNP, (Pathogen_SNP &  Host_SNP) | (!Pathogen_SNP & !Host_SNP))
    Ct_value_data_SNP_Both_None$Pathogen_Host <- Ct_value_data_SNP_Both_None$Pathogen_SNP & Ct_value_data_SNP_Both_None$Host_SNP
    
    print(summary(lm('Ct_value_int ~ HIV_status + symptoms_duration_cough_duration + TB_RF_smoking + Host_SNP*Pathogen_SNP', 
                     data = Ct_value_data_SNP_Both_None)))
    print(summary(lm('Ct_value_int ~ HIV_status + symptoms_duration_cough_duration + TB_RF_smoking + Host_SNP', 
                     data = Ct_value_data_SNP)))
    print(summary(lm('Ct_value_int ~ HIV_status + symptoms_duration_cough_duration  + TB_RF_smoking + Pathogen_SNP', 
                     data = Ct_value_data_SNP)))

  }else if(analysis == 'Host'){
    host_genotype <- snpStats::read.plink(bed = '/mnt/data3/zmxu/G2G_TB/data/Genotyping_WGS/TBDAR.WGS.Imputed.GWASReady',select.snps = host_snp)
    PCs <- data.table::fread('/mnt/data3/zmxu/G2G_TB/data/Genotyping_WGS/TBDAR.WGS.Imputed.GWASReady.eigenvec') %>% dplyr::select(IID = V1,PC1 = V3,PC2 = V4,PC3 = V5)
    PCs$PATIENT_ID <- sapply(PCs$IID,function(x) strsplit(x=x,split = '_')[[1]][2])
    PCs$PATIENT_ID <- as.numeric(sapply(PCs$PATIENT_ID,function(x) gsub(x=x,pattern = 'Fellay.',replacement = '')))
    
    
    host_dosage <- as(host_genotype$genotypes, Class = 'numeric')
    #Flip allele, REF allele as 0
    host_dosage <- abs(host_dosage-2)
    host_dosage <- as.data.frame(host_dosage)
    host_dosage$PATIENT_ID <- sapply(host_genotype$fam$pedigree,function(x) strsplit(x=x,split = '_')[[1]][2])
    host_dosage$PATIENT_ID <- as.numeric(sapply(host_dosage$PATIENT_ID,function(x) gsub(x=x,pattern = 'Fellay.',replacement = '')))
    
    metadata <- data.table::fread('../../data/pheno/metadata_Sinergia_final_all_pats_corrected.txt',stringsAsFactors = F) %>% dplyr::select(PATIENT_ID,TB_score,Ct_value,Xray_score,age,
                                                                                                                       HIV_status,patient_sex,symptoms_duration_cough_duration,TB_RF_smoking) %>% dplyr::distinct(PATIENT_ID,.keep_all = T)
    metadata$HIV_status[metadata$HIV_status == ''] <- NA
    metadata$HIV_status[metadata$HIV_status == 'infected'] <- 'positive'
    
    metadata$TB_RF_smoking[metadata$TB_RF_smoking == 'unknown'] <- NA
    metadata$TB_RF_smoking[metadata$TB_RF_smoking == ''] <- NA
    
    metadata$Xray_score_binary <- 0
    #High severity with X-ray > 70
    metadata$Xray_score_binary[metadata$Xray_score > 70] <- 1
    metadata$Xray_score_binary[is.na(metadata$Xray_score)] <- NA
    
    metadata$Ct_value_log <- log10(as.numeric(metadata$Ct_value))
    
    #metadata <- metadata %>% dplyr::select(-Ct_value,-Xray_score)
    
    jned_data <- dplyr::left_join(host_dosage,metadata) %>% dplyr::left_join(PCs)
    data.table::fwrite(jned_data,'../../data/Genotyping_WGS/metadata_tbdarfinal_humanavaliable.txt',sep = '\t',quote = F,na = 'NA',row.names = F)

    #TB Score Regression
    print(coeftest(polr(glue::glue("factor(TB_score) ~ HIV_status + symptoms_duration_cough_duration + TB_RF_smoking +  patient_sex + {colnames(host_dosage)[1]}"), 
                        data = jned_data,Hess = T)))
    print(nrow(na.omit(jned_data %>% dplyr::select(TB_score,HIV_status,symptoms_duration_cough_duration,TB_RF_smoking,patient_sex,one_of(colnames(host_dosage)[1])))))
    
    print(coeftest(polr(glue::glue("factor(TB_score) ~ HIV_status + symptoms_duration_cough_duration + TB_RF_smoking +  patient_sex + {colnames(host_dosage)[2]}"), 
                        data = jned_data,Hess = T)))
    print(nrow(na.omit(jned_data %>% dplyr::select(TB_score,HIV_status,symptoms_duration_cough_duration,TB_RF_smoking,patient_sex,one_of(colnames(host_dosage)[2])))))
    
    
    #X-ray score 
    print(summary(glm(glue::glue("Xray_score_binary  ~  age + {colnames(host_dosage)[1]}"), 
                      data = jned_data)))
    print(nrow(na.omit(jned_data %>% dplyr::select(Xray_score_binary,age,one_of(colnames(host_dosage)[1])))))
    
    print(summary(glm(glue::glue("Xray_score_binary  ~  age + {colnames(host_dosage)[2]}"), 
                                 data = jned_data)))
    print(nrow(na.omit(jned_data %>% dplyr::select(Xray_score_binary,age,one_of(colnames(host_dosage)[2])))))
    
    #Ct value
    print(summary(lm(glue::glue("Ct_value_log ~ HIV_status + symptoms_duration_cough_duration  + {colnames(host_dosage)[1]}"), 
                     data = jned_data)))
    print(nrow(na.omit(jned_data %>% dplyr::select(Ct_value_log,HIV_status,symptoms_duration_cough_duration,one_of((colnames(host_dosage)[1]))))))
          
    print(summary(lm(glue::glue("Ct_value_log ~ HIV_status + symptoms_duration_cough_duration  + {colnames(host_dosage)[2]}"), 
                     data = jned_data)))
    print(nrow(na.omit(jned_data %>% dplyr::select(Ct_value_log,HIV_status,symptoms_duration_cough_duration,one_of((colnames(host_dosage)[2]))))))
    
  }else if(analysis == 'Pathogen'){
    metadata <- data.table::fread('../../data/Mtb/metadata_tbdarfinal_bacavailable.txt')
    
    metadata$HIV_status[metadata$HIV_status == ''] <- NA
    metadata$HIV_status[metadata$HIV_status == 'infected'] <- 'positive'
    
    metadata$TB_RF_smoking[metadata$TB_RF_smoking == 'unknown'] <- NA
    metadata$TB_RF_smoking[metadata$TB_RF_smoking == ''] <- NA
    
    metadata$Xray_score_binary <- 0
    #High severity with X-ray > 70
    metadata$Xray_score_binary[metadata$Xray_score > 70] <- 1
    metadata$Xray_score_binary[is.na(metadata$Xray_score)] <- NA
    
    metadata$Ct_value_log <- log10(as.numeric(metadata$Ct_value))
    
    metadata <- metadata %>% dplyr::select(-Ct_value,-Xray_score)
    
    #TB Score Regression
    print(coeftest(polr(glue::glue("factor(TB_score) ~ HIV_status + symptoms_duration_cough_duration + TB_RF_smoking +  patient_sex + Mutation_3388671"), 
                        data = metadata,Hess = T)))
    print(nrow(na.omit(metadata %>% dplyr::select(TB_score,HIV_status,symptoms_duration_cough_duration,TB_RF_smoking,patient_sex,Mutation_3388671))))
    
    print(coeftest(polr(glue::glue("factor(TB_score) ~ HIV_status + symptoms_duration_cough_duration + TB_RF_smoking +  patient_sex + Mutation_2626678"), 
                        data = metadata,Hess = T)))
    print(nrow(na.omit(metadata %>% dplyr::select(TB_score,HIV_status,symptoms_duration_cough_duration,TB_RF_smoking,patient_sex,Mutation_2626678))))
    
    #X-ray score 
    print(summary(glm(glue::glue("Xray_score_binary  ~  age + Mutation_3388671"), 
                      data = metadata)))
    print(nrow(na.omit(metadata %>% dplyr::select(Xray_score_binary,age,Mutation_3388671))))
    print(summary(glm(glue::glue("Xray_score_binary  ~  age + Mutation_2626678"), 
                      data = metadata)))
    print(nrow(na.omit(metadata %>% dplyr::select(Xray_score_binary,age,Mutation_2626678))))
    
    #Ct value
    print(summary(lm(glue::glue("Ct_value_log ~ HIV_status + symptoms_duration_cough_duration  + Mutation_3388671"), 
                     data = metadata)))
    print(nrow(na.omit(metadata %>% dplyr::select(Ct_value_log,HIV_status,symptoms_duration_cough_duration,Mutation_3388671))))
    
    print(summary(lm(glue::glue("Ct_value_log ~ HIV_status + symptoms_duration_cough_duration  + Mutation_2626678"), 
                     data = metadata)))
    print(nrow(na.omit(metadata %>% dplyr::select(Ct_value_log,HIV_status,symptoms_duration_cough_duration,Mutation_2626678))))
    
    
  }
  
  
  # ggplot2::ggplot(data = pathogen_data_jned %>% dplyr::filter(!is.na(Pathogen_SNP)),aes(x=factor(Pathogen_SNP),y=TB_score)) + geom_boxplot() + stat_n_text() +   stat_compare_means(method = "anova") 
  # ggplot2::ggplot(data = pathogen_data_jned %>% dplyr::filter(!is.na(Pathogen_SNP)),aes(x=factor(Pathogen_SNP),y=Xray_score)) + geom_boxplot() + stat_n_text() +   stat_compare_means(method = "anova") 
  # ggplot2::ggplot(data = pathogen_data_jned %>% dplyr::filter(!is.na(Pathogen_SNP)),aes(x=factor(Pathogen_SNP),y=Ct_value)) + geom_boxplot() + stat_n_text() +   stat_compare_means(method = "anova") 
  # 
  # ggplot2::ggplot(data = pathogen_data_filt,aes_string(x='Group',y=pheno)) + geom_boxplot() + stat_n_text() + stat_compare_means(method = "anova") 
  
}
# CheckPheno(G2G_Obj = readRDS('../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds'),
#            host_snp = 'rs12151990',
#            AA_variant = 'Rv2348c_Rv2348c:2626678:p.Ile101Met',
#            lineage = 'ALL',
#            pheno = 'TB_score')

# CheckPheno(G2G_Obj = readRDS('../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds'),
#            host_snp = 'rs75769176',
#            AA_variant = 'fixA_Rv3029c:3388671:p.Thr67Met',
#            lineage = 'ALL',
#            pheno = 'TB_score',analysis = 'Both')

# CheckPheno(G2G_Obj = readRDS('../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds'),
#            host_snp = c('rs12151990','rs75769176'),
#            pheno = 'TB_score',analysis = 'Host')

# CheckPheno(G2G_Obj = readRDS('../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds'),
#            AA_variant = c('fixA_Rv3029c:3388671:p.Thr67Met','Rv2348c_Rv2348c:2626678:p.Ile101Met'),
#            pheno = 'TB_score',analysis = 'Pathogen')
