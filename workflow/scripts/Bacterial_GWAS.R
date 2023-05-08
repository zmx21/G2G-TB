Mtb_AAs = data.table::fread('/mnt/data1/yildiran/TBDAR/Mtb/Mtb_AAs.txt') %>% dplyr::select(-FID,-IID)
Metadata <- data.table::fread('/mnt/data1/yildiran/TBDAR/Metadata/TBDAR_Metadata_Complete.tsv')
Metadata$Xray_score_binary <- 0
Metadata$Xray_score_binary[Metadata$Xray_score > 71] <- 1
Metadata$Ct_value <- log10(Metadata$Ct_value)

TB_P <- rep(NA,ncol(Mtb_AAs))
Xray_P <-rep(NA,ncol(Mtb_AAs))
Ct_value_P <- rep(NA,ncol(Mtb_AAs))

for(i in 1:ncol(Mtb_AAs)){
  cur_bac_var <- Mtb_AAs[,..i]
  colnames(cur_bac_var) <- 'Bac_SNP'
  cur_bac_var_name <- colnames(Mtb_AAs)[i]
  
  TB_score_data <- na.omit(cbind(cur_bac_var,dplyr::select(Metadata,TB_score,age,HIV_status,patient_sex,symptoms_duration_cough_duration,TB_RF_smoking)))
  Xray_score_data <- na.omit(cbind(cur_bac_var,dplyr::select(Metadata,Xray_score_binary,age,HIV_status,patient_sex,symptoms_duration_cough_duration,TB_RF_smoking)))
  Ct_value_data <- na.omit(cbind(cur_bac_var,dplyr::select(Metadata,Ct_value,age,HIV_status,patient_sex,symptoms_duration_cough_duration,TB_RF_smoking)))

  TB_mdl <- summary(lm("TB_score ~ HIV_status + symptoms_duration_cough_duration + TB_RF_smoking +  patient_sex + Bac_SNP" ,data = TB_score_data))
  ct_mdl <- summary(lm(glue::glue("Ct_value ~ HIV_status + symptoms_duration_cough_duration  + Bac_SNP") , data = Ct_value_data))
  Xray_mdl <- summary(glm(glue::glue("Xray_score_binary ~  Bac_SNP") , data = Xray_score_data))
 
  TB_P[i] <-  TB_mdl$coefficients['Bac_SNP','Pr(>|t|)']
  Ct_value_P[i] <-  ct_mdl$coefficients['Bac_SNP','Pr(>|t|)']
  Xray_P[i] <-  Xray_mdl$coefficients['Bac_SNP','Pr(>|t|)']
  
}
