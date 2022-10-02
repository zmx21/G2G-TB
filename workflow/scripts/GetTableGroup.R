library(snpStats)
GetAASNPData <- function(AA_Hit,SNP_Hit,G2G_Obj){
  Mtb_AA <- G2G_Obj$aa_matrix_full[,AA_Hit]
  Mtb_AA[Mtb_AA==2] <- 1
  
  host_genotype <- snpStats::read.plink(bed = paste0(G2G_Obj$host_path['ALL'],'.bed'),select.snps = SNP_Hit)
  host_dosage <- as(host_genotype$genotypes, Class = 'numeric')
  host_dosage <- host_dosage[,1]
  host_dosage[host_dosage==2] <- paste0(host_genotype$map$allele.2,host_genotype$map$allele.2)
  host_dosage[host_dosage==1] <- paste0(host_genotype$map$allele.1,host_genotype$map$allele.2)
  host_dosage[host_dosage==0] <- paste0(host_genotype$map$allele.1,host_genotype$map$allele.1)
  host_dosage_df <- data.frame(FAM_ID = names(host_dosage),Host_SNP = host_dosage)
  colnames(host_dosage_df) <- c('FAM_ID',SNP_Hit)
  
  df <- data.frame(G_NUMBER = rownames(G2G_Obj$aa_matrix_full),Mtb_AA = Mtb_AA)
  colnames(df) <- c('G_NUMBER',AA_Hit)
  
  df <- df %>% dplyr::left_join(G2G_Obj$both_IDs_to_keep$ALL,by=c('G_NUMBER'='G_NUMBER')) %>% dplyr::left_join(host_dosage_df,by = c('FAM_ID'='FAM_ID')) %>% dplyr::rename(IID=FAM_ID) %>%
    dplyr::relocate(G_NUMBER,PATIENT_ID,IID,LINEAGE)
  return(df)
}
CalculateHostDistance <- function(PC_Df,from_centroid = F){
  if(!from_centroid){
    dist_matrix <- matrix(NA,nrow = nrow(PC_Df),ncol = nrow(PC_Df))
    rownames(dist_matrix) <- PC_Df$IID
    colnames(dist_matrix) <- PC_Df$IID
    
    for(i in 1:nrow(dist_matrix)){
      for(j in 1:ncol(dist_matrix)){
        dist_matrix[i,j] <- sqrt((PC_Df$PC1[i] - PC_Df$PC1[j])^2 + (PC_Df$PC2[i] - PC_Df$PC2[j])^2)
      }
    }
    return(dist_matrix)
  }else if(from_centroid) {
    dist_df <- data.frame(IID=character(),Host_Dist=numeric())
    centroid_x <- mean(PC_Df$PC1)
    centroid_y <- mean(PC_Df$PC2)
    
    for(i in 1:nrow(PC_Df)){
      dist_df <- rbind(dist_df,data.frame(IID=PC_Df$IID[i],Host_Dist = sqrt((PC_Df$PC1[i] - centroid_x)^2 + (PC_Df$PC2[i] - centroid_y)^2)))
    }
    return(dist_df)
  }
}
PC_Dist_Centroid <- CalculateHostDistance(readRDS('../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds')$host_PCs$ALL,from_centroid = T)
PC_Dist_Centroid$Host_Dist <- PC_Dist_Centroid$Host_Dist / sd(PC_Dist_Centroid$Host_Dist)

hiv_status <- data.table::fread('../data/pheno/metadata_Sinergia_final_dataset_human_bac_genome_available_QCed.txt') %>% dplyr::select(PATIENT_ID,HIV_status)
Rv2348c_df <- GetAASNPData("Rv2348c_Rv2348c:2626678:p.Ile101Met",'rs12151990',readRDS('../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds')) %>% dplyr::left_join(hiv_status) %>% dplyr::left_join(PC_Dist_Centroid)
data.table::fwrite(Rv2348c_df,'../experiment_grouping/rs12151990_Rv2348c-Rv2348c-2626678-p.Ile101Met.tsv',col.names = T,row.names = F,quote = F,na = 'NA',sep = '\t')
fixA_df <- GetAASNPData("fixA_Rv3029c:3388671:p.Thr67Met",'rs75769176',readRDS('../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds')) %>% dplyr::left_join(hiv_status) %>% dplyr::left_join(PC_Dist_Centroid)
data.table::fwrite(fixA_df,'../experiment_grouping/rs75769176_fixA-Rv3029c-3388671-p.Thr67Met.tsv',col.names = T,row.names = F,quote = F,na = 'NA',sep = '\t')

# PC_Dist <- CalculateHostDistance(readRDS('../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds')$host_PCs$ALL)
# write.table(PC_Dist,file = '../experiment_grouping/Host_PC_Dist.txt')

Bac_Dist <- data.table::fread('../data/Mtb/tbdarfinal_bacavailable.distance_matrix')
bac_row_names <- Bac_Dist$V1
Bac_Dist <- Bac_Dist[,-'V1']
Bac_Dist <- as.matrix(Bac_Dist)
rownames(Bac_Dist) <- bac_row_names

# Approach 1
Group1_Approach1_Rv2348 <- dplyr::filter(Rv2348c_df,`Rv2348c_Rv2348c:2626678:p.Ile101Met`==1 & rs12151990=='TC' & LINEAGE == 'L4') %>% dplyr::arrange(Host_Dist) %>% dplyr::select(-HIV_status,-IID,-LINEAGE) %>% 
  dplyr::rename(Rv2348c_I101M=`Rv2348c_Rv2348c:2626678:p.Ile101Met`) %>% dplyr::arrange(Host_Dist)
Group1_Approach1_fixA <- dplyr::filter(fixA_df,`fixA_Rv3029c:3388671:p.Thr67Met`==1 & rs75769176=='GA' & LINEAGE == 'L3') %>% dplyr::arrange(Host_Dist) %>% dplyr::select(-HIV_status,-IID,-LINEAGE) %>% 
  dplyr::rename(fixA_T67M=`fixA_Rv3029c:3388671:p.Thr67Met`) %>% dplyr::arrange(Host_Dist)

Group2_Approach1_Rv2348  <- dplyr::filter(Rv2348c_df,`Rv2348c_Rv2348c:2626678:p.Ile101Met`==0 & rs12151990=='TC') %>% dplyr::arrange(Host_Dist) %>% dplyr::filter(HIV_status == 'negative') %>% dplyr::select(-HIV_status,-IID) %>% 
  dplyr::filter(G_NUMBER %in% rownames(Bac_Dist) & LINEAGE == 'L4') %>% dplyr::rename(Rv2348c_I101M=`Rv2348c_Rv2348c:2626678:p.Ile101Met`) %>% dplyr::select(-LINEAGE)
Group2_Approach1_Rv2348$Bac_Dist <- sapply(Group2_Approach1_Rv2348$G_NUMBER,function(x) mean(Bac_Dist[x,Group1_Approach1_Rv2348$G_NUMBER]))
Group2_Approach1_Rv2348 %>% dplyr::arrange(Bac_Dist,Host_Dist)
Group2_Approach1_fixA  <- dplyr::filter(fixA_df,`fixA_Rv3029c:3388671:p.Thr67Met`==0 & rs75769176=='GA') %>% dplyr::arrange(Host_Dist) %>% dplyr::filter(HIV_status == 'negative') %>% dplyr::select(-HIV_status,-IID) %>% 
  dplyr::filter(G_NUMBER %in% rownames(Bac_Dist) & LINEAGE == 'L3') %>% dplyr::rename(fixA_T67M=`fixA_Rv3029c:3388671:p.Thr67Met`) %>% dplyr::select(-LINEAGE)
Group2_Approach1_fixA$Bac_Dist <- sapply(Group2_Approach1_fixA$G_NUMBER,function(x) mean(Bac_Dist[x,Group1_Approach1_fixA$G_NUMBER]))
Group2_Approach1_fixA %>% dplyr::arrange(Bac_Dist,Host_Dist)

Group3_Approach1_Rv2348  <- dplyr::filter(Rv2348c_df,`Rv2348c_Rv2348c:2626678:p.Ile101Met`==1 & rs12151990=='CC') %>% dplyr::arrange(Host_Dist) %>% dplyr::filter(HIV_status == 'negative') %>% dplyr::select(-HIV_status,-IID) %>% 
  dplyr::filter(G_NUMBER %in% rownames(Bac_Dist) & LINEAGE == 'L4') %>% dplyr::rename(Rv2348c_I101M=`Rv2348c_Rv2348c:2626678:p.Ile101Met`) %>% dplyr::select(-LINEAGE)
Group3_Approach1_Rv2348$Bac_Dist <- sapply(Group3_Approach1_Rv2348$G_NUMBER,function(x) mean(Bac_Dist[x,Group1_Approach1_Rv2348$G_NUMBER]))
Group3_Approach1_Rv2348 %>% dplyr::arrange(Bac_Dist,Host_Dist)
Group3_Approach1_fixA  <- dplyr::filter(fixA_df,`fixA_Rv3029c:3388671:p.Thr67Met`==1 & rs75769176=='AA') %>% dplyr::arrange(Host_Dist) %>% dplyr::filter(HIV_status == 'negative') %>% dplyr::select(-HIV_status,-IID) %>% 
  dplyr::filter(G_NUMBER %in% rownames(Bac_Dist) & LINEAGE == 'L3') %>% dplyr::rename(fixA_T67M=`fixA_Rv3029c:3388671:p.Thr67Met`) %>% dplyr::select(-LINEAGE)
Group3_Approach1_fixA$Bac_Dist <- sapply(Group3_Approach1_fixA$G_NUMBER,function(x) mean(Bac_Dist[x,Group1_Approach1_fixA$G_NUMBER]))
Group3_Approach1_fixA %>% dplyr::arrange(Bac_Dist,Host_Dist)

Group4_Approach1_Rv2348  <- dplyr::filter(Rv2348c_df,`Rv2348c_Rv2348c:2626678:p.Ile101Met`==0 & rs12151990=='CC') %>% dplyr::arrange(Host_Dist) %>% dplyr::filter(HIV_status == 'negative') %>% dplyr::select(-HIV_status,-IID) %>% 
  dplyr::filter(G_NUMBER %in% rownames(Bac_Dist) & LINEAGE == 'L4') %>% dplyr::rename(Rv2348c_I101M=`Rv2348c_Rv2348c:2626678:p.Ile101Met`) %>% dplyr::select(-LINEAGE)
Group4_Approach1_Rv2348$Bac_Dist <- sapply(Group4_Approach1_Rv2348$G_NUMBER,function(x) mean(Bac_Dist[x,Group1_Approach1_Rv2348$G_NUMBER]))
Group4_Approach1_Rv2348<- Group4_Approach1_Rv2348 %>% dplyr::arrange(Bac_Dist,Host_Dist)
Group4_Approach1_Rv2348[1:50,]
Group4_Approach1_fixA  <- dplyr::filter(fixA_df,`fixA_Rv3029c:3388671:p.Thr67Met`==0 & rs75769176=='AA') %>% dplyr::arrange(Host_Dist) %>% dplyr::filter(HIV_status == 'negative') %>% dplyr::select(-HIV_status,-IID) %>% 
  dplyr::filter(G_NUMBER %in% rownames(Bac_Dist) & LINEAGE == 'L3') %>% dplyr::rename(fixA_T67M=`fixA_Rv3029c:3388671:p.Thr67Met`) %>% dplyr::select(-LINEAGE)
Group4_Approach1_fixA$Bac_Dist <- sapply(Group4_Approach1_fixA$G_NUMBER,function(x) mean(Bac_Dist[x,Group1_Approach1_fixA$G_NUMBER]))
Group4_Approach1_fixA <- Group4_Approach1_fixA %>% dplyr::arrange(Bac_Dist,Host_Dist)
Group4_Approach1_fixA[1:50,]

# Approach 2

Group1_Approach2_Rv2348 <- Group1_Approach1_Rv2348
Group4_Approach2_Rv2348 <- dplyr::filter(Rv2348c_df,`Rv2348c_Rv2348c:2626678:p.Ile101Met`==0 & rs12151990=='CC') %>% dplyr::arrange(Host_Dist) %>% dplyr::filter(HIV_status == 'negative') %>% dplyr::select(-HIV_status,-IID) %>% 
  dplyr::filter(G_NUMBER %in% rownames(Bac_Dist) & LINEAGE == 'L4') %>% dplyr::rename(Rv2348c_I101M=`Rv2348c_Rv2348c:2626678:p.Ile101Met`) %>% dplyr::select(-LINEAGE)
Group4_Approach2_Rv2348$Bac_Dist <- sapply(Group4_Approach2_Rv2348$G_NUMBER,function(x) mean(Bac_Dist[x,Group1_Approach2_Rv2348$G_NUMBER]))
Group4_Approach2_Rv2348 %>% dplyr::arrange(Bac_Dist,Host_Dist)



