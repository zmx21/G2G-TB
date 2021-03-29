GetParsedMapping <- function(DATA_DIR = '~/G2G_TB/data/'){
  #Get genotyped patient IDs
  fam_file_genotyped <- data.table::fread(glue::glue('{DATA_DIR}Genotyping/TB_DAR_Genotyping_WGS_Merged/WGS.genotyped.merged.fam'))
  genotyped_IDs <- sapply(fam_file_genotyped$V2,function(x) ifelse(grepl(x=x,pattern = '-'),strsplit(x=x,split = '-')[[1]][1],x))
  genotyped_IDs <- sapply(genotyped_IDs,function(x) strsplit(x=x,split = "_")[[1]][2])
  genotyped_IDs <- sapply(genotyped_IDs,function(x) ifelse(grepl(x=x,pattern = 'Fellay\\.'),strsplit(x=x,split = 'Fellay\\.')[[1]][2],x))
  
  
  #Old meta data file, sent on 27/28/2020
  meta1 <- data.table::fread(glue::glue("{DATA_DIR}/pheno/metadata_bloodbacgenomes_27082020.txt")) %>% dplyr::filter(PATIENT_ID %in% genotyped_IDs)
  #New meta data file, sent on 16/10/2020
  meta2 <- readxl::read_xlsx(glue::glue("{DATA_DIR}/pheno/metadata_2020-10-16_17-36-52.xlsx")) %>% dplyr::filter(`PATIENT ID` %in% genotyped_IDs)
  #Third meta data file, sent on 12/02/2021
  meta3 <- data.table::fread(glue::glue("{DATA_DIR}/pheno/metadata_combined_genomes_022021.txt")) %>% dplyr::filter(PATIENT_ID %in% genotyped_IDs)
  
  #IDs based on VCF files 
  VCF_DIR <- '~/G2G_TB/data/Mtb/Mtb_VCF_files/full/'
  VCF_Files <- dir(VCF_DIR)[grepl(pattern = 'var.homo',x=dir(VCF_DIR))]
  VCF_IDs <- sapply(VCF_Files,function(x) strsplit(x=x,split = '\\.')[[1]][1])
  
  #Patient ID to G Number, based on metadata1
  meta1_mapping <- dplyr::select(meta1,PATIENT_ID,G_NUMBER,LINEAGE.1=`G_NUMBER/LINEAGE`)
  meta1_mapping$PATIENT_ID <- as.character(meta1_mapping$PATIENT_ID)
  
  #Patient ID to G Number, based on metadata2
  meta2_mapping <- dplyr::select(meta2,PATIENT_ID=`PATIENT ID`,G_NUMBER=`G NUMBER`)
  
  #Patient ID to G Number, based on metadata3
  meta3_mapping <- dplyr::select(meta3,PATIENT_ID,G_NUMBER,LINEAGE.3=Lineage)
  meta3_mapping$PATIENT_ID <- as.character(meta3_mapping$PATIENT_ID)
  
  
  
  #Discrepancy between the two mapping
  jned_mapping <- dplyr::left_join(meta2_mapping,meta1_mapping,by=c('PATIENT_ID'='PATIENT_ID')) %>% dplyr::left_join(meta3_mapping,by=c('PATIENT_ID'='PATIENT_ID')) %>% dplyr::rename(G_NUMBER.2 = G_NUMBER.x,
                                                                                                                                                                                      G_NUMBER.1 = G_NUMBER.y,
                                                                                                                                                                                      G_NUMBER.3 = G_NUMBER)
  disrep_mapping <- jned_mapping[which(jned_mapping$G_NUMBER.1 != jned_mapping$G_NUMBER.2),]
  #nrow(disrep_mapping) #85 patients with disrepancy in mapping
  #View(disrep_mapping)
  
  #Check which mapping VCF files correspond to
  #length(intersect(VCF_IDs,disrep_mapping$G_NUMBER.1)) # Most VCFs (85) correspond to ID in meta1
  #length(intersect(VCF_IDs,disrep_mapping$G_NUMBER.2)) #Some VCFs (5) correspond to ID in meta2
  
  #disrep_mapping[(disrep_mapping$G_NUMBER.1 %in% VCF_IDs & disrep_mapping$G_NUMBER.2 %in% VCF_IDs),] #For some patients (5 in total), both G Number correspond to a VCF file
  #disrep_mapping[(disrep_mapping$G_NUMBER.1 %in% VCF_IDs & !disrep_mapping$G_NUMBER.2 %in% VCF_IDs),] #For some patients (5 in total), both G Number correspond to a VCF file
  
  #Check which mapping tree file correspond to
  tree <- ape::read.tree(glue::glue("{DATA_DIR}/pheno/RAxML_bestTree.Sinergia_MTB_1088_052020_rerooted"))
  tree_tips <- tree$tip.label
  length(intersect(tree_tips,disrep_mapping$G_NUMBER.1)) # Some VCFs (5) correspond to ID in meta1
  length(intersect(tree_tips,disrep_mapping$G_NUMBER.2)) #Most VCFs (82) correspond to ID in meta2
  #disrep_mapping[(disrep_mapping$G_NUMBER.1 %in% tree_tips & disrep_mapping$G_NUMBER.1 %in% tree_tips),] #For the same patients as the VCF files (5 in total), both G Number correspond to a tree tip. 
  
  parsed_mapping <- rbind(dplyr::filter(jned_mapping,G_NUMBER.3 %in% VCF_IDs) %>% dplyr::select(PATIENT_ID,G_NUMBER = G_NUMBER.3,LINEAGE = LINEAGE.3),
                          dplyr::filter(jned_mapping,!G_NUMBER.3 %in% VCF_IDs & G_NUMBER.2 %in% VCF_IDs) %>% dplyr::select(PATIENT_ID,G_NUMBER = G_NUMBER.2,LINEAGE = LINEAGE.1),
                          dplyr::filter(jned_mapping,!G_NUMBER.3 %in% VCF_IDs & !G_NUMBER.2 %in% VCF_IDs & G_NUMBER.1 %in% VCF_IDs) %>% dplyr::select(PATIENT_ID,G_NUMBER = G_NUMBER.1,LINEAGE = LINEAGE.1))
  return(parsed_mapping)
}
parsed_mapping <- GetParsedMapping()
