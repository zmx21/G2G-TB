tb_dar_pheno <- data.table::fread('../../pheno/metadata_Sinergia_final_all_pats_corrected.txt') %>% dplyr::mutate(PATIENT_ID = as.character(PATIENT_ID))

WGS_Sex_Check <- data.table::fread('../../WGS/WGS_Host_Data/plink.sexcheck') %>% dplyr::select(PLINK_ID=IID,SNPSEX) %>% dplyr::mutate(Source = 'WGS')
WGS_Sex_Check$SNPSEX <- ifelse(WGS_Sex_Check$SNPSEX == 1, 'MALE','FEMALE')
WGS_Sex_Check$PATIENT_ID <- sapply(WGS_Sex_Check$PLINK_ID,function(x) gsub(x=x,pattern = 'WGS_Fellay.',replacement = ''))

Batch1_Sex_Check <- data.table::fread('../../Genotyping/TB_DAR_Genotyping/Batch1/Sex_Est_Table.txt') %>% dplyr::select(PATIENT_ID=`Sample ID`,SNPSEX = Gender) %>% dplyr::mutate(SNPSEX = toupper(SNPSEX),Source = 'Batch1')
Batch1_Sex_Check$PATIENT_ID <- sapply(Batch1_Sex_Check$PATIENT_ID,function(x) strsplit(x=x,split = '-')[[1]][1])
Batch1_Sex_Check$PATIENT_ID <- sapply(Batch1_Sex_Check$PATIENT_ID,function(x) gsub(x=x,pattern = ' V01 B1',replacement = ''))
Batch1_Sex_Check$PLINK_ID <- paste0('Batch1_',Batch1_Sex_Check$PATIENT_ID)

Batch2_Sex_Check <- data.table::fread('../../Genotyping/TB_DAR_Genotyping/Batch2/Sex_Est_Table.txt') %>% dplyr::select(PATIENT_ID=`Sample ID`,SNPSEX = Gender) %>% dplyr::mutate(SNPSEX = toupper(SNPSEX),Source = 'Batch2')
Batch2_Sex_Check$PATIENT_ID <- sapply(Batch2_Sex_Check$PATIENT_ID,function(x) strsplit(x=x,split = '-')[[1]][1])
Batch2_Sex_Check$PATIENT_ID <- sapply(Batch2_Sex_Check$PATIENT_ID,function(x) gsub(x=x,pattern = 'TZ ',replacement = ''))
Batch2_Sex_Check$PLINK_ID <- paste0('Batch2_',Batch2_Sex_Check$PATIENT_ID)

sex_comparison <- dplyr::left_join(rbind(WGS_Sex_Check,Batch1_Sex_Check,Batch2_Sex_Check),tb_dar_pheno %>% dplyr::select(patient_sex,PATIENT_ID) %>% dplyr::mutate(patient_sex = toupper(patient_sex)),by = c('PATIENT_ID'='PATIENT_ID'))

sex_mismatch <- dplyr::filter(sex_comparison,SNPSEX != patient_sex & SNPSEX != 'UNKNOWN')

data.table::fwrite(sex_mismatch %>% dplyr::select(FID=PLINK_ID,IID=PLINK_ID),'Sex_Mismatch.txt',sep = ' ',row.names = F,quote = F)
