library(dplyr)
GetParsedMapping <- function(Genotyping_DIR,Pheno_DIR,Phylo_Tree_Path){
  #Get genotyped patient IDs
  fam_file_genotyped <- data.table::fread(glue::glue('{Genotyping_DIR}TB_DAR_Genotyping_WGS_Merged/WGS.genotyped.merged.fam'))
  genotyped_IDs <- sapply(fam_file_genotyped$V2,function(x) ifelse(grepl(x=x,pattern = '-'),strsplit(x=x,split = '-')[[1]][1],x))
  genotyped_IDs <- sapply(genotyped_IDs,function(x) strsplit(x=x,split = "_")[[1]][2])
  genotyped_IDs <- sapply(genotyped_IDs,function(x) ifelse(grepl(x=x,pattern = 'Fellay\\.'),strsplit(x=x,split = 'Fellay\\.')[[1]][2],x))
  
  
  #Old meta data file, sent on 27/28/2020
  meta1 <- data.table::fread(glue::glue("{Pheno_DIR}metadata_bloodbacgenomes_27082020.txt")) %>% dplyr::filter(PATIENT_ID %in% genotyped_IDs)
  #New meta data file, sent on 16/10/2020
  meta2 <- readxl::read_xlsx(glue::glue("{Pheno_DIR}metadata_2020-10-16_17-36-52.xlsx")) %>% dplyr::filter(`PATIENT ID` %in% genotyped_IDs)
  #Third meta data file, sent on 12/02/2021
  meta3 <- data.table::fread(glue::glue("{Pheno_DIR}metadata_combined_genomes_022021_inclXray.txt")) %>% dplyr::filter(PATIENT_ID %in% genotyped_IDs)
  
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
  tree <- ape::read.tree(Phylo_Tree_Path)
  tree_tips <- tree$tip.label
  length(intersect(tree_tips,disrep_mapping$G_NUMBER.1)) # Some VCFs (5) correspond to ID in meta1
  length(intersect(tree_tips,disrep_mapping$G_NUMBER.2)) #Most VCFs (82) correspond to ID in meta2
  #disrep_mapping[(disrep_mapping$G_NUMBER.1 %in% tree_tips & disrep_mapping$G_NUMBER.1 %in% tree_tips),] #For the same patients as the VCF files (5 in total), both G Number correspond to a tree tip. 
  
  parsed_mapping <- rbind(dplyr::filter(jned_mapping,G_NUMBER.3 %in% VCF_IDs) %>% dplyr::select(PATIENT_ID,G_NUMBER = G_NUMBER.3,LINEAGE = LINEAGE.3),
                          dplyr::filter(jned_mapping,!G_NUMBER.3 %in% VCF_IDs & G_NUMBER.2 %in% VCF_IDs) %>% dplyr::select(PATIENT_ID,G_NUMBER = G_NUMBER.2,LINEAGE = LINEAGE.1),
                          dplyr::filter(jned_mapping,!G_NUMBER.3 %in% VCF_IDs & !G_NUMBER.2 %in% VCF_IDs & G_NUMBER.1 %in% VCF_IDs) %>% dplyr::select(PATIENT_ID,G_NUMBER = G_NUMBER.1,LINEAGE = LINEAGE.1))
  return(parsed_mapping)
}

GetMAC <- function(AA_Matrix_filt){
  MAC <- apply(AA_Matrix_filt,2,function(x) {
    no_na <- x[!is.na(x)]
    if(length(no_na) == 0){
      return(0)
    }
    no_na_binary <- ifelse(no_na == 0,0,1)
    return(ifelse(length(unique(no_na_binary)) != 1,min(table(no_na_binary)),length(no_na_binary) - min(table(no_na_binary))))
  })
  return(MAC)
}
ConstructpPCA <- function(AA_Table,Phylo_Tree_Path,ID_Mapping_df,pPCA = T,mac_thresh = 5){
  tree <- ape::read.tree(Phylo_Tree_Path)
  ID_Mapping_df_filt <- dplyr::filter(ID_Mapping_df,G_NUMBER %in% tree$tip.label & !is.na(LINEAGE))
  
  #Separate based on LINEAGES
  unique_lineages <- sort(unique(ID_Mapping_df_filt$LINEAGE))
  unique_lineages <- c(unique_lineages,unlist(sapply(2:(length(unique_lineages)-1),function(x) apply(combn(unique_lineages,x),2,function(y) paste0(sort(y),collapse = '_')))))
  lineage_samples <- lapply(unique_lineages,function(x) dplyr::filter(ID_Mapping_df_filt,LINEAGE %in% strsplit(x,split = '_')[[1]])$G_NUMBER)
  #Add in index which cover all LINEAGES
  lineage_samples <- c(list(ID_Mapping_df_filt$G_NUMBER),lineage_samples)
  names(lineage_samples) <- c('ALL',unique_lineages)
  
  vir_pca <- vector(mode = 'list',length = length(lineage_samples))
  names(vir_pca) <- c('ALL',unique_lineages)
  for(i in 1:length(lineage_samples)){
    tree_subset <-ape::keep.tip(tree, lineage_samples[[i]])
    AA_Table_filt <- AA_Table[tree_subset$tip.label,]
    MAC <- GetMAC(AA_Table_filt)
    
    if(pPCA){
      vir_4d <- phylobase::phylo4d(tree_subset, AA_Table_filt[,MAC > mac_thresh])
      vir_pca[[i]] <- adephylo::ppca(
        vir_4d,
        scale = FALSE,
        center = T,
        scannf = F,
        nfposi = 10,
        method = "Abouheif"
      )
    }else{
      vir_pca[[i]] <- prcomp(AA_Table_filt[,MAC > mac_thresh],center = T,scale = F)
    }
  }
  return(vir_pca)
}
#Filter AA variants according to MAC and Group variants according to lineage
#If a variant appears in multiple lineages, do unstratified analysis.
#If a variant appears in only 1 lineage, do stratified analysis.
FilterAAMatrix <- function(AA_Matrix,Lineage_Df,MAC_Thresh){
  #Keep samples where Host data exists
  AA_Matrix_filt <- AA_Matrix[Lineage_Df$ALL$G_NUMBER,,drop=FALSE]
  #Keep variants which pass unstratified MAC threshold
  MAC_AA_Matrix_filt <- GetMAC(AA_Matrix_filt)
  AA_Matrix_filt <- AA_Matrix_filt[,MAC_AA_Matrix_filt > MAC_Thresh,drop=FALSE]
  #Remove Missing Variants/Genes
  AA_Matrix_filt <- AA_Matrix_filt[,!apply(AA_Matrix_filt,2,function(x) all(is.na(x))),drop=FALSE]
  
  #Assign variants to lineages
  strat_table <- lapply(1:ncol(AA_Matrix_filt),function(i) table(ifelse(AA_Matrix_filt[,i]==0,0,1),Lineage_Df$ALL$LINEAGE))
  strat_assng <- vector(mode = 'character',length = length(strat_table))
  for(i in 1:length(strat_table)){
    cur_strat <- strat_table[[i]]
    MAC_by_lineage <- apply(cur_strat,2,function(x) min(x,na.rm = T))
    #If a variant is perfectly stratified by lineage and there are no variability within a lineage, exclude from analysis
    if(sum(MAC_by_lineage > MAC_Thresh) == 0){
      strat_assng[i] <- NA
      #If a variant only show variability in a single lineage, assign it to that lineage
    }else if (sum(MAC_by_lineage > MAC_Thresh) == 1){
      strat_assng[i] <- names(MAC_by_lineage)[which(MAC_by_lineage > MAC_Thresh)]
    }else{
      strat_assng[i] <- paste0(sort(names(MAC_by_lineage)[which(MAC_by_lineage > MAC_Thresh)]),collapse = '_')
      if(strat_assng[i] == 'L1_L2_L3_L4'){
        strat_assng[i] <- 'ALL'
      }
    }
  }
  
  #Construct AA Matrix for each lineage
  uniq_assgn <- sort(unique(strat_assng))
  uniq_assgn <- uniq_assgn[!is.na(uniq_assgn)]
  
  #Exclude ALL lineage (Will be run seperately anyways)
  #uniq_assgn <- uniq_assgn[uniq_assgn != paste0(sort(unique(Lineage_Df$ALL$LINEAGE)),collapse = '_')] 
  #For each gene, assign to lineage groups
  AA_Matrix_filt_by_lineage <- lapply(uniq_assgn,function(x) AA_Matrix_filt[Lineage_Df[[x]]$G_NUMBER,which(strat_assng == x),drop=FALSE])
  names(AA_Matrix_filt_by_lineage) <- uniq_assgn
  
  return(AA_Matrix_filt_by_lineage)
}


SetUpG2G <- function(OUT_DIR,Genotyping_DIR,Pheno_DIR,Var_Tbl_Path,Phylo_Tree_Path,Mtb_Nuc_Out,Ref_Panel,Host_MAF,Pathogen_MAC_pPCA,Pathogen_MAC_AA_Lineage,
                     Pathogen_MAC_AA,VCF_Path,n_cores){

  #Get mapping file with consensus samples
  parsed_mapping <- GetParsedMapping(Genotyping_DIR,Pheno_DIR,Phylo_Tree_Path)
  
  #Get nucleotide variant matrix
  nuc_matrix <- data.table::fread(Mtb_Nuc_Out)
  colnames(nuc_matrix) <- sapply(colnames(nuc_matrix),function(x) strsplit(x=x,split = ']')[[1]][2])
  nuc_matrix_trans <- nuc_matrix[,-c('CHROM','POS','REF','ALT')]
  nuc_matrix_trans <- t(as.matrix(nuc_matrix_trans))
  nuc_matrix_trans[nuc_matrix_trans==2] <- 1
  rownames(nuc_matrix_trans) <- colnames(nuc_matrix[,-c('CHROM','POS','REF','ALT')])
  colnames(nuc_matrix_trans) <- paste0(nuc_matrix$POS,'_',nuc_matrix$REF,'_',nuc_matrix$ALT)
  
  #Calculate Mtb pPCs
  vir_pPCA <- ConstructpPCA(nuc_matrix_trans,Phylo_Tree_Path,parsed_mapping,Pathogen_MAC_pPCA)
  vir_pPCs <- lapply(vir_pPCA,function(x) cbind(G_NUMBER=rownames(x$li),as.data.frame(x$li)) %>% dplyr::left_join(parsed_mapping))
  
  #Get consensus samples 
  host_IDs <- data.table::fread(glue::glue("{VCF_Path}/TB_DAR_Imputed_GWAS.fam"))$V2
  WGS_index <- sapply(host_IDs,function(x) grepl(x=x,pattern = 'WGS_Fellay'))
  if(any(WGS_index)){
    host_genotyped_IDs <- host_IDs[!WGS_index]
    host_genotyped_IDs_simple <- sapply(host_genotyped_IDs,function(x) strsplit(x=x,split = '_')[[1]][2])
    host_genotyped_IDs_simple <- sapply(host_genotyped_IDs_simple,function(x) ifelse(grepl(x=x,pattern = '-'),strsplit(x=x,split = '-')[[1]][1],x))
    
    host_WGS_IDs_simple <- sapply(host_IDs[WGS_index],function(x) strsplit(x=x,split = 'WGS_Fellay\\.')[[1]][2])
    
    host_IDs_simple <- rep(NA,length(host_IDs))
    host_IDs_simple[!WGS_index] <- host_genotyped_IDs_simple
    host_IDs_simple[WGS_index] <- host_WGS_IDs_simple
    names(host_IDs_simple) <- host_IDs
    
  }else{
    host_genotyped_IDs <- sapply(host_IDs,function(x) strsplit(x=x,split = '_')[[1]][2])
    host_IDs_simple <- sapply(host_genotyped_IDs,function(x) ifelse(grepl(x=x,pattern = '-'),strsplit(x=x,split = '-')[[1]][1],x))
  }
  
  #Initialize variables to store
  both_IDs_to_keep <- vector(mode = 'list',length = length(vir_pPCs)); names(both_IDs_to_keep) <- names(vir_pPCs)
  host_PCs <- vector(mode = 'list',length = length(vir_pPCs)); names(host_PCs) <- names(vir_pPCs)
  covars <- vector(mode = 'list',length = length(vir_pPCs)); names(covars) <- names(vir_pPCs)
  
  #Get Host PCs
  host_PCs_raw <- data.table::fread(glue::glue("{VCF_Path}/TB_DAR_Imputed_GWAS_PCA.eigenvec"))
  
  #Get Covars
  covars_discrete <- data.table::fread(glue::glue("{VCF_Path}/covars_discrete")) 
  covars_numeric <- data.table::fread(glue::glue("{VCF_Path}/covars_numeric")) 
  covars_raw <- dplyr::inner_join(covars_discrete,covars_numeric)
  
  #Loop through lineages, store PLINK files, covars, and host PCs
  for(i in 1:length(vir_pPCs)){
    system(glue::glue('mkdir -p {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/'))
    pathogen_IDs_to_keep <- dplyr::filter(parsed_mapping,G_NUMBER %in% vir_pPCs[[i]]$G_NUMBER)
    both_IDs_to_keep[[i]] <- dplyr::filter(pathogen_IDs_to_keep,PATIENT_ID %in% host_IDs_simple)
    both_IDs_to_keep[[i]]$FAM_ID <- host_IDs[match(both_IDs_to_keep[[i]]$PATIENT_ID,host_IDs_simple)]
    write(both_IDs_to_keep[[i]]$FAM_ID,file = glue::glue("{OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt"))
    
    #Write out PLINK files with samples to keep
    system(
      glue::glue(
        "plink2 --bfile {VCF_Path}/TB_DAR_Imputed_GWAS --keep {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --indiv-sort f {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --make-bed --out {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_Imputed_G2G_Raw"
      )
    )
    
    #Apply MAF Filter
    system(
      glue::glue(
        "plink --bfile {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_Imputed_G2G_Raw --impute-sex --make-bed --out {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_Imputed_G2G_Raw2"
      )
    )
    
    #Impute Sex
    system(
      glue::glue(
        "plink2 --bfile {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_Imputed_G2G_Raw2 --maf {Host_MAF} --indiv-sort f {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --make-bed --out {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_Imputed_G2G"
      )
    )
    
    #Write out HLA PLINK files with samples to keep
    system(
      glue::glue(
        "plink2 --bfile {VCF_Path}/TB_DAR_HLA_Alleles --keep {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --indiv-sort f {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --make-bed --out {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_HLA_Alleles_G2G"
      )
    )
    system(
      glue::glue(
        "plink2 --bfile {VCF_Path}/TB_DAR_HLA_AA --keep {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --indiv-sort f {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --make-bed --out {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_HLA_AA_G2G"
      )
    )
    
    system(glue::glue("rm {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_Imputed_G2G_Raw*"))
    #Order and filter host PCs
    host_PCs[[i]] <- dplyr::left_join(data.frame('IID' = both_IDs_to_keep[[i]]$FAM_ID),host_PCs_raw[,-1],by=c('IID'='V2'))
    colnames(host_PCs[[i]]) <- c('IID',paste0('PC',seq(1,ncol(host_PCs[[i]])-1,1)))
    
    #Order and filter pPCs
    vir_pPCs[[i]] <- dplyr::left_join(data.frame(G_NUMBER=both_IDs_to_keep[[i]]$G_NUMBER),vir_pPCs[[i]])
    
    #Order and filter covars
    covars[[i]] <- dplyr::left_join(data.frame(IID=both_IDs_to_keep[[i]]$FAM_ID),covars_raw) %>% dplyr::select(-contains("PC"))
  }
  #Unstratified AA matrix
  aa_matrix <- readRDS(Var_Tbl_Path) 
  aa_matrix_raw <- aa_matrix[both_IDs_to_keep[['ALL']]$G_NUMBER,,drop=FALSE]
  aa_matrix_raw <- aa_matrix_raw[,GetMAC(aa_matrix_raw) > Pathogen_MAC_AA,drop=FALSE]
  
  
  #Filter AA Matrix, decide for each variant whether to do stratified or stratified analysis
  aa_matrix_filt <- FilterAAMatrix(aa_matrix_raw,both_IDs_to_keep,MAC_Thresh = Pathogen_MAC_AA_Lineage) 
  
  #Path to VCF file (for each lineage)
  host_path <- glue::glue("{OUT_DIR}/{Ref_Panel}/LINEAGE_{c(names(vir_pPCs),'ALL')}/TB_DAR_Imputed_G2G")
  names(host_path) <- c(names(vir_pPCs),'ALL')
  
  #Load TB Score
  tb_score <- data.table::fread(glue::glue("{VCF_Path}/tb_score"))
  tb_score_by_lineage <- lapply(both_IDs_to_keep,function(x) dplyr::left_join(x %>% dplyr::select(IID=FAM_ID),tb_score,c('IID'='IID')) %>% dplyr::relocate(FID,IID,tb_score))
  
  #Load X-Ray Score
  xray_score <- data.table::fread(glue::glue("{VCF_Path}/xray_score"))
  xray_score_by_lineage <- lapply(both_IDs_to_keep,function(x) dplyr::left_join(x %>% dplyr::select(IID=FAM_ID),xray_score,c('IID'='IID')) %>% dplyr::relocate(FID,IID,xray_score))
  
  return(list(host_PCs=host_PCs,
              vir_pPCs=vir_pPCs,
              covars = covars,
              tb_score=tb_score_by_lineage,
              xray_score = xray_score_by_lineage,
              aa_matrix_filt=aa_matrix_filt,
              aa_matrix_raw=aa_matrix_raw,
              nuc_matrix = nuc_matrix_trans,
              both_IDs_to_keep=both_IDs_to_keep,
              host_path = host_path,
              raw_pPCA = vir_pPCA))
}


args <- commandArgs(trailingOnly = TRUE) 
OUT_DIR <- args[[1]]
Genotyping_DIR <- args[[2]]
Pheno_DIR <- args[[3]]
Var_Tbl_Path <- args[[4]]
Phylo_Tree_Path <- args[[5]]
Mtb_Nuc_Out <- args[[6]]
Ref_Panel <- args[[7]]
Host_MAF <- as.numeric(args[[8]])
Pathogen_MAC_pPCA <- as.numeric(args[[9]])
Pathogen_MAC_AA_Lineage <- as.numeric(args[[10]])
Pathogen_MAC_AA <- as.numeric(args[[11]])
VCF_Path <- args[[12]]
n_cores <- as.numeric(args[[13]])

G2G_Obj <- SetUpG2G(OUT_DIR,Genotyping_DIR,Pheno_DIR,Var_Tbl_Path,Phylo_Tree_Path,Mtb_Nuc_Out,Ref_Panel,Host_MAF,Pathogen_MAC_pPCA,Pathogen_MAC_AA_Lineage,Pathogen_MAC_AA,VCF_Path,n_cores)
saveRDS(G2G_Obj,glue::glue("{OUT_DIR}/{Ref_Panel}/G2G_Obj.rds"))

Pheno_DIR <- '~/G2G_TB/data/pheno/'
meta_final <- data.table::fread('~/G2G_TB/data/pheno/metadata_final.txt')
VCF_DIR <- '~/G2G_TB/data/Mtb/Mtb_VCF_files/full/'
VCF_Files <- dir(VCF_DIR)[grepl(pattern = 'var.homo',x=dir(VCF_DIR))]
VCF_IDs <- sapply(VCF_Files,function(x) strsplit(x=x,split = '\\.')[[1]][1])

meta_final_mapping <- dplyr::select(meta_final,PATIENT_ID,G_NUMBER)
meta_final_mapping$PATIENT_ID <- as.character(meta_final_mapping$PATIENT_ID)

WGS_IDs_QC <- system('~/Software/bcftools query -l ~/G2G_TB/data/WGS/WGS_Host_Data/WGS_Fellay.hg38.joint.118h-1947437863.genotyped.renamed.consensus.filt.vcf.gz',intern = T)
WGS_IDs_QC <- sapply(WGS_IDs_QC,function(x) strsplit(x=x,split = 'WGS_Fellay\\.')[[1]][2])

genotyped_IDs <- data.table::fread('~/G2G_TB/data/Genotyping/TB_DAR_Genotyping/archive/PLINK_Raw/Fellay_Plates1to12.samples',header = F)
genotyped_IDs <- genotyped_IDs$V1
genotyped_IDs <- sapply(genotyped_IDs,function(x) strsplit(x=x,split = '-')[[1]][1])

genotyped_IDs_QC <- system('~/Software/bcftools query -l ~/G2G_TB/data/Genotyping/TB_DAR_Genotyping/PLINK_QCed_Autosomes_Top/VCF/Fellay_0620.GRCh37.vcf.gz',intern = T)
genotyped_IDs_QC <- sapply(genotyped_IDs_QC,function(x) strsplit(x,split = '_')[[1]][2])
genotyped_IDs_QC <- sapply(genotyped_IDs_QC,function(x) strsplit(x,split = '-')[[1]][1])

batch_3 <- data.table::fread('~/G2G_TB/data/Genotyping/TB_DAR_Genotyping/Batch2/Fellay_0122_corrected/Samples Table.txt')
genotyped_IDs_batch_3 <- batch_3$`Sample ID`
genotyped_IDs_batch_3 <- gsub(genotyped_IDs_batch_3,pattern = 'TZ ',replacement = '')
genotyped_IDs_batch_3 <- sapply(genotyped_IDs_batch_3,function(x) strsplit(x=x,split = '-')[[1]][1])

all_Ids <- data.frame(PATIENT_ID = unique(c(genotyped_IDs,genotyped_IDs_batch_3,genotyped_IDs_QC,WGS_IDs_QC)))


status_df <- all_Ids %>%  dplyr::left_join(meta_final_mapping)

status_df$Batch_1 <- ifelse(status_df$PATIENT_ID %in% WGS_IDs_QC,'PASS',NA)
status_df$Batch_2 <- ifelse(status_df$PATIENT_ID %in% genotyped_IDs,ifelse(status_df$PATIENT_ID %in% genotyped_IDs_QC,'PASS','FAIL'),NA)
status_df$Batch_3 <- ifelse(status_df$PATIENT_ID %in% genotyped_IDs_batch_3,'PASS',NA)
status_df$TB_Genome <- ifelse(status_df$G_NUMBER %in% VCF_IDs,T,F)
data.table::fwrite(status_df,file = '~/G2G_TB/data/Sinergia_Final_Host_Bac_Status.txt',sep = '\t',col.names = T,quote = F,na = 'NA')

