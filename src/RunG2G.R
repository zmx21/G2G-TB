library(pbmcapply)
library(dplyr)
library(parallel)
library(mltools)
library(data.table)
source('~/TB_GWAS/src/TB_GWAS.R')
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

ConstructpPCA <- function(AA_Table,Phylo_Tree_Path,ID_Mapping_df,pPCA = T,mac_thresh = 5){
  tree <- ape::read.tree(Phylo_Tree_Path)
  ID_Mapping_df_filt <- dplyr::filter(ID_Mapping_df,G_NUMBER %in% tree$tip.label & !is.na(LINEAGE))
  
  #Separate based on LINEAGES
  unique_lineages <- sort(unique(ID_Mapping_df_filt$LINEAGE))
  lineage_samples <- lapply(unique_lineages,function(x) dplyr::filter(ID_Mapping_df_filt,LINEAGE == x)$G_NUMBER)
  #Add in index which cover all LINEAGES
  lineage_samples <- c(list(ID_Mapping_df_filt$G_NUMBER),lineage_samples)
  names(lineage_samples) <- c('ALL',unique_lineages)
  
  vir_pca <- vector(mode = 'list',length = length(lineage_samples))
  names(vir_pca) <- c('ALL',unique_lineages)
  for(i in 1:length(lineage_samples)){
    tree_subset <-ape::keep.tip(tree, lineage_samples[[i]])
    AA_Table_filt <- AA_Table[tree_subset$tip.label,]
    MAC <- apply(AA_Table_filt,2,function(x) ifelse(length(unique(x)) != 1,min(table(x)),nrow(AA_Table_filt) - min(table(x))))
    
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
  AA_Matrix_filt <- AA_Matrix[Lineage_Df$ALL$G_NUMBER,]
  #Keep variants which pass unstratified MAC threshold
  MAC_AA_Matrix_filt <-  apply(AA_Matrix_filt,2,function(x) ifelse(length(unique(x)) != 1,min(table(x)),nrow(AA_Matrix_filt) - min(table(x))))
  AA_Matrix_filt <- AA_Matrix_filt[,MAC_AA_Matrix_filt > MAC_Thresh]
  
  #Assign variants to lineages
  strat_table <- lapply(1:ncol(AA_Matrix_filt),function(i) table(AA_Matrix_filt[,i],Lineage_Df$ALL$LINEAGE))
  strat_assng <- vector(mode = 'character',length = length(strat_table))
  for(i in 1:length(strat_table)){
    cur_strat <- strat_table[[i]]
    MAC_by_lineage <- apply(cur_strat,2,function(x) min(x))
    #If a variant is perfectly stratified by lineage and there are no variability within a lineage, exclude from analysis
    if(sum(MAC_by_lineage > MAC_Thresh) == 0){
      strat_assng[i] <- NA
      #If a variant only show variability in a single lineage, assign it to that lineage
    }else if (sum(MAC_by_lineage > MAC_Thresh) == 1){
      strat_assng[i] <- names(MAC_by_lineage)[which(MAC_by_lineage > MAC_Thresh)]
    }else{
      strat_assng[i] <- 'ALL'
    }
  }
  
  #Construct AA Matrix for each lineage
  uniq_assgn <- sort(unique(strat_assng))
  uniq_assgn <- uniq_assgn[!is.na(uniq_assgn)]
  AA_Matrix_filt_by_lineage <- lapply(uniq_assgn,function(x) AA_Matrix_filt[Lineage_Df[[x]]$G_NUMBER,which(strat_assng == x)])
  names(AA_Matrix_filt_by_lineage) <- uniq_assgn
  
  return(AA_Matrix_filt_by_lineage)
}


SetUpG2G <- function(DATA_DIR,SOFTWARE_DIR,OUT_DIR,Ref_Panel = 'AFGR',Host_MAF,Pathogen_MAC_pPCA=5,Pathogen_MAC_AA=10,VCF_Path = NULL,SIFT = F){
  PLINK = '/home/zmxu/Software/plink'
  
  Mtb_Data <- paste0(DATA_DIR,"Mtb/Mtb_VCF_files/full/")
  Mtb_Out <- paste0(OUT_DIR,"Mtb_AA_Tbl/")
  Mtb_Out_Syn <- paste0(OUT_DIR,"Mtb_Syn_Tbl/")
  
  Phylo_Tree_Path <- '~/G2G_TB/data/Mtb/fasttree_TBDar_022021_1239'
  Mtb_Nuc_Out <- paste0(DATA_DIR,"Mtb/Mtb_VCF_files/merged/merged.var.homo.SNPs.dosage")
  
  #Get mapping file with consensus samples
  parsed_mapping <- GetParsedMapping(DATA_DIR)
  
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
  
  #Generate QCed file based on specified reference panel
  if(is.null(VCF_Path)){
    VCF_Path <- glue::glue("{DATA_DIR}Genotyping/TB_DAR_{Ref_Panel}_Imputed/{Ref_Panel}.imputed.vcf.gz")
  }
  # SetUpGWAS(VCF_Path,glue::glue("{OUT_DIR}/{Ref_Panel}/"),maf_thresh=Host_MAF,covars_discrete_to_incl = c('Patient_Sex','HIV_Status','TB_RF_Smoking','BMI','Patient_Household_Size'),covars_numeric_to_incl = c('Age'),col_names=T)
  
  #Get consensus samples 
  host_IDs <- data.table::fread(glue::glue("{OUT_DIR}/{Ref_Panel}/TB_DAR_Imputed_GWAS.fam"))$V2
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
  host_PCs_raw <- data.table::fread(glue::glue("{OUT_DIR}/{Ref_Panel}/TB_DAR_Imputed_GWAS_PCA.eigenvec"))
  
  #Get Covars
  covars_discrete <- data.table::fread(glue::glue("{OUT_DIR}/{Ref_Panel}/covars_discrete")) 
  covars_numeric <- data.table::fread(glue::glue("{OUT_DIR}/{Ref_Panel}/covars_numeric")) 
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
        "{PLINK}2 --bfile {OUT_DIR}/{Ref_Panel}/TB_DAR_Imputed_GWAS --keep {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --indiv-sort f {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --make-bed --out {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_Imputed_G2G_Raw"
      )
    )
    
    #Apply MAF Filter
    system(
      glue::glue(
        "{PLINK} --bfile {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_Imputed_G2G_Raw --impute-sex --make-bed --out {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_Imputed_G2G_Raw2"
      )
    )
    
    #Impute Sex
    system(
      glue::glue(
        "{PLINK}2 --bfile {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_Imputed_G2G_Raw2 --maf {Host_MAF} --indiv-sort f {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --make-bed --out {OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_Imputed_G2G"
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
  
  #Filter AA Matrix, decide for each variant whether to do stratified or stratified analysis
  aa_matrix <- readRDS(glue::glue("{Mtb_Out}AA_Matrix.rds")) 
  aa_matrix_filt <- FilterAAMatrix(aa_matrix,both_IDs_to_keep,MAC_Thresh = Pathogen_MAC_AA) 
  
  #Unstratified AA matrix
  aa_matrix_raw <- aa_matrix[both_IDs_to_keep[['ALL']]$G_NUMBER,]
  aa_matrix_MAC <- apply(aa_matrix_raw,2,function(x) ifelse(length(unique(x)) != 1,min(table(x)),nrow(aa_matrix_raw) - min(table(x))))
  aa_matrix_raw <- aa_matrix_raw[,aa_matrix_MAC > Pathogen_MAC_AA]

  #Gene Burden Matrix
  if(SIFT){
    gene_burden <- readRDS(glue::glue("{Mtb_Out_Syn}Gene_Burden_SIFT.rds")) 
  }else{
    gene_burden <- readRDS(glue::glue("{Mtb_Out_Syn}Gene_Burden.rds")) 
  }
  gene_burden_non_syn <- gene_burden$non_syn_burden[both_IDs_to_keep[['ALL']]$G_NUMBER,]
  gene_burden_syn <- gene_burden$syn_burden[both_IDs_to_keep[['ALL']]$G_NUMBER,]
  
  #Path to VCF file (for each lineage)
  host_path <- glue::glue("{OUT_DIR}/{Ref_Panel}/LINEAGE_{names(vir_pPCs)}/TB_DAR_Imputed_G2G")
  names(host_path) <- names(vir_pPCs)
  
  #Load TB Score
  tb_score <- data.table::fread(glue::glue("{OUT_DIR}/{Ref_Panel}/tb_score"))
  tb_score_by_lineage <- lapply(both_IDs_to_keep,function(x) dplyr::left_join(x %>% dplyr::select(IID=FAM_ID),tb_score,c('IID'='IID')) %>% dplyr::relocate(FID,IID,tb_score))
  
  return(list(host_PCs=host_PCs,
              vir_pPCs=vir_pPCs,
              covars = covars,
              tb_score=tb_score_by_lineage,
              aa_matrix_filt=aa_matrix_filt,
              aa_matrix_raw=aa_matrix_raw,
              nuc_matrix = nuc_matrix_trans,
              both_IDs_to_keep=both_IDs_to_keep,
              host_path = host_path,
              raw_pPCA = vir_pPCA,
              gene_burden_non_syn=gene_burden_non_syn,
              gene_burden_syn = gene_burden_syn))
}

GetResults <- function(OUT_DIR,suffix = 'glm.logistic.hybrid',p_thresh=5e-8,n_cores=5,is_interaction = F,is_ordinal = F,tool = 'PLINK'){
  all_files <- dir(OUT_DIR)
  result_files <- all_files[sapply(all_files,function(x) grepl(pattern = suffix,x=x))]
  if(!is_interaction){
    if((tool == 'PLINK' | tool == 'PLINK-FIRTH') & !grepl(x=suffix,pattern = 'gz')){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0("awk '{ if (NR == 1 || $13 <= ",p_thresh,") {print} }' ",OUT_DIR,x)),mc.cores = n_cores)
    }else if ((tool == 'PLINK' | tool == 'PLINK-FIRTH') & grepl(x=suffix,pattern = 'gz')){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0('zcat ',OUT_DIR,x," | awk '{ if (NR == 1 || $13 <= ",p_thresh,") {print} }' ")),mc.cores = n_cores)
    }else if(tool == 'GMMAT-SCORE'){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0('zcat ',OUT_DIR,x," | awk '{ if (NR == 1 || $11 <= ",p_thresh,") {print} }' ")),mc.cores = n_cores)
    }
  }else{
    if(is_ordinal){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0("zcat ",OUT_DIR,x," | awk -F ',' '{ if ($6 <= ",p_thresh,") {print} }' ")),mc.cores = n_cores)
    }else{
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0("zcat ",OUT_DIR,x," | grep 'ADDx' | awk '{ if ($12 <= ",p_thresh,") {print} }' ")),mc.cores = n_cores)
    }
  }
  names(results) <- result_files
  return(results)
}

RunG2G <- function(G2G_Obj,SOFTWARE_DIR,OUT_DIR,Ref_Panel,tool = 'GCTA',n_PC = 5,n_pPC = 6,n_cores = 20,covars_to_incl = c(),lineage = c(),debug=F,chunks = c(),chr = c(),test_snp = c(),model = NA,X_Chr=F){
  OUT_PATH <- glue::glue("{OUT_DIR}/{Ref_Panel}")
  system(glue::glue("mkdir -p {OUT_PATH}"))
  
  GCTA = '/home/zmxu/Software/gcta'
  PLINK = '/home/zmxu/Software/plink'
  bcftools = '/home/zmxu/Software/bcftools'
  regenie = '/home/zmxu/Software/regenie_v2.0.1.gz_x86_64_Linux'
  OUT_PATH <- glue::glue("{OUT_PATH}/{tool}/")
  system(glue::glue("mkdir -p {OUT_PATH}"))
  if(length(covars_to_incl) > 0){
    OUT_PATH <- glue::glue("{OUT_PATH}/PC_{n_PC}_pPC_{n_pPC}_cov_{paste0(sapply(covars_to_incl,function(x) gsub(x=x,pattern='_',replacement='')),collapse = '-')}/")
  }else{
    OUT_PATH <- glue::glue("{OUT_PATH}/PC_{n_PC}_pPC_{n_pPC}/")
  }
  system(glue::glue("mkdir -p {OUT_PATH}"))
  
  if (tool == 'PLINK' | tool == 'PLINK-FIRTH' | tool == 'SAIGE' | tool == 'GMMAT-SCORE' | tool == 'GMMAT-WALD'){
    if(length(lineage)==0){
      lineages_to_run <- unique(names(G2G_Obj$aa_matrix_filt))
    }else{
      lineages_to_run <- lineage
    }
    for(i in 1:length(lineages_to_run)){
      cur_lineage <- lineages_to_run[i]
      OUT_PATH_Lineage <- glue::glue("{OUT_PATH}/LINEAGE_{cur_lineage}/")
      system(glue::glue("mkdir -p {OUT_PATH_Lineage}"))
      system(glue::glue("mkdir -p {OUT_PATH_Lineage}/tmp"))
      
      #Get IDD and FID, ensure they're in correct order
      fam_file <- data.table::fread(glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"))
      colnames(fam_file)[1:2] <- c('FID','IID')
      if(!all(fam_file$IID == G2G_Obj$both_IDs_to_keep[[cur_lineage]]$FAM_ID)){
        stop('Sample order incorrect')
      }
      
      #Write out AA matrix
      cur_aa_matrix_filt <- G2G_Obj$aa_matrix_filt[[cur_lineage]]
      cur_aa_Matrix_mat <- cbind(fam_file[,c(1,2)],as.matrix(cur_aa_matrix_filt))
      colnames(cur_aa_Matrix_mat) <- c('FID','IID',colnames(cur_aa_matrix_filt))
      
      if(debug){
        cur_aa_Matrix_mat <- cur_aa_Matrix_mat[,1:5]
        data.table::fwrite(cur_aa_Matrix_mat,col.names = T,row.names = F,sep = ' ',file = glue::glue("{OUT_PATH_Lineage}/tmp/AA_outcome.txt"),na = 'NA',quote = F)
      }else{
        data.table::fwrite(cur_aa_Matrix_mat,col.names = T,row.names = F,sep = ' ',file = glue::glue("{OUT_PATH_Lineage}/tmp/AA_outcome.txt"),na = 'NA',quote = F)
      }
      
      
      #Store PCs and covars
      host_PCs <- G2G_Obj$host_PCs[[cur_lineage]]
      host_PCs <- dplyr::select(host_PCs,'IID',paste0('PC',1:n_PC))
      
      if(n_pPC != 0){
        pPCs <- G2G_Obj$vir_pPCs[[cur_lineage]]
        pPCs <-dplyr::select(pPCs,IID,paste0('PC',1:n_pPC))
        colnames(pPCs) <- c('IID',paste0('pPC',1:n_pPC))
      }
      if(length(covars_to_incl) > 0){
        covars_num_discrete <- G2G_Obj$covars[[cur_lineage]]
        covars_num_discrete <- dplyr::select(covars_num_discrete,'IID',covars_to_incl)
        covars_num_discrete[is.na(covars_num_discrete)] <- 'NONE'
        
        if(n_pPC != 0){
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(pPCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(covars_num_discrete,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }else{
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(covars_num_discrete,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }
      }else{
        if(n_pPC != 0){
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(pPCs,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }else{
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }
      }
      
      #Change categorial binary covariates to numeric (0 and 1)
      cat_covars <- which(sapply(3:ncol(covars),function(q) any(is.character(as.vector(t(covars[,..q]))))))
      if(length(cat_covars) > 0){
        for(col in (cat_covars+2)){
          covars[,col] <- as.integer(as.factor(t(covars[,..col]))) - 1
        }
      }
      
      data.table::fwrite(covars,col.names = F,row.names = F,sep = ' ',file = glue::glue("{OUT_PATH_Lineage}/tmp/plink-covars.txt"),na = 'NA',quote = F)
      
      #Run association study for each AA variant in the current lineage
      AA_Matrix_No_ID <- cur_aa_Matrix_mat[,-c(1,2)]
      #Run association study for each AA variant in the current lineage
      if(length(chunks) == 0){
        chunks <- 1:ncol(AA_Matrix_No_ID)
      }
      if(max(chunks) > ncol(AA_Matrix_No_ID)){
        chunks <- chunks[1]:ncol(AA_Matrix_No_ID)
      }
      
      #Preprocess for GMATT
      if((tool == 'GMMAT-WALD' | tool == 'GMMAT-SCORE') & !file.exists(glue::glue("{OUT_PATH_Lineage}/tmp/output/GEMMA.sXX.txt"))){
        #Run GEMMA to make GRM (https://www.xzlab.org/software/GEMMAmanual.pdf)
        cur_dir <- getwd()
        setwd(glue::glue("{OUT_PATH_Lineage}/tmp/"))
        fam_file <- data.table::fread(glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"),header = F)
        fam_file$V6 <- 0 #GEMMA dislikes -9 label (but doesn't use pheno info)
        data.table::fwrite(fam_file,glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"),col.names = F,sep = ' ',append = F)
        
        system(glue::glue("{SOFTWARE_DIR}/gemma-0.98.1-linux-static -bfile {G2G_Obj$host_path[cur_lineage]} -gk 2 -o /GEMMA"))
        setwd(cur_dir)
      }
      
      for(k in chunks){
        cur_pathogen_variant <- colnames(AA_Matrix_No_ID)[k]
        if(tool == 'PLINK'){
          if(is.na(model)){
            if(X_Chr){
              tryCatch(system(
                glue::glue(
                  "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --glm cc-residualize hide-covar --covar-variance-standardize --1 --pheno {OUT_PATH_Lineage}/tmp/AA_outcome.txt --pheno-col-nums {k+2} --covar {OUT_PATH_Lineage}/tmp/plink-covars.txt --out {OUT_PATH_Lineage}{cur_pathogen_variant}"
                )
              ))
              
            }else{
              tryCatch(system(
                glue::glue(
                  "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --no-sex --glm cc-residualize hide-covar --covar-variance-standardize --1 --pheno {OUT_PATH_Lineage}/tmp/AA_outcome.txt --pheno-col-nums {k+2} --covar {OUT_PATH_Lineage}/tmp/plink-covars.txt --out {OUT_PATH_Lineage}{cur_pathogen_variant}"
                )
              ))
              
            }
            
            system(glue::glue('pigz --fast {OUT_PATH_Lineage}{cur_pathogen_variant}*.hybrid'))
          }else{
            tryCatch(system(
              glue::glue(
                "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --no-sex --glm {model} cc-residualize hide-covar --covar-variance-standardize --1 --pheno {OUT_PATH_Lineage}/tmp/AA_outcome.txt --pheno-col-nums {k+2} --covar {OUT_PATH_Lineage}/tmp/plink-covars.txt --out {OUT_PATH_Lineage}{cur_pathogen_variant}.{model}"
              )
            ))
            
          }
        }else if(tool == 'PLINK-FIRTH'){
          tryCatch(system(
            glue::glue(
              "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --no-sex --glm firth hide-covar --1 --pheno {OUT_PATH_Lineage}/tmp/AA_outcome.txt --pheno-col-nums {k+2} --covar {OUT_PATH_Lineage}/tmp/plink-covars.txt --out {OUT_PATH_Lineage}{cur_pathogen_variant}"
            )
          ))
        }else if(tool == 'SAIGE'){
          saige_pheno <- dplyr::left_join(data.frame(cur_aa_Matrix_mat),covars,by=c('IID'='IID','FID'='FID')) %>% dplyr::select(-FID) %>% dplyr::rename(ID = IID)
          data.table::fwrite(saige_pheno,col.names = T,row.names = F,quote = F,na = 'NA',sep = ' ',glue::glue("{OUT_PATH_Lineage}/tmp/SAIGE-pheno.txt"))
          nullGLMM_stat <- tryCatch(expr = {
            SAIGE::fitNULLGLMM(plinkFile = G2G_Obj$host_path[cur_lineage],
                               phenoFile = glue::glue("{OUT_PATH_Lineage}/tmp/SAIGE-pheno.txt"),
                               phenoCol = cur_pathogen_variant,
                               traitType = 'binary',
                               invNormalize = F,
                               covarColList = setdiff(colnames(covars),c("IID","FID")),
                               sampleIDColinphenoFile = "ID",
                               LOCO = F,
                               nThreads = n_cores,
                               outputPrefix = glue::glue('{OUT_PATH_Lineage}{cur_pathogen_variant}'))
            T
          },error = function(e){
            write("NULL MODEL ERROR",glue::glue('{OUT_PATH_Lineage}{cur_pathogen_variant}.error'),append = T)
            write(as.character(e),glue::glue('{OUT_PATH_Lineage}{cur_pathogen_variant}.error'),append = T)
            return(F)},
          warning = function(w) {
            write("NULL MODEL WARNING",glue::glue('{OUT_PATH_Lineage}{cur_pathogen_variant}.warning'),append = T)
            write(as.character(w),glue::glue('{OUT_PATH_Lineage}{cur_pathogen_variant}.warning'),append = T)
            return(F)})
          
          system(
            glue::glue(
              "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --export vcf-4.2 bgz --out {G2G_Obj$host_path[cur_lineage]}"
            )
          )
          
          system(
            glue::glue(
              "{bcftools} index -t {G2G_Obj$host_path[cur_lineage]}.vcf.gz"
            )
          )
          
          system(
            glue::glue(
              "{bcftools} query -l {G2G_Obj$host_path[cur_lineage]}.vcf.gz > {G2G_Obj$host_path[cur_lineage]}.vcf.gz.sample"
            )
          )
          #Run through all chr if no specified. 
          if(length(chr) == 0){
            chr <- 1:22
          }
          err_counter <- mclapply(chr,function(j){
            SPA_test_stat <- tryCatch(expr = {
              SAIGE::SPAGMMATtest(vcfFile = glue::glue("{G2G_Obj$host_path[cur_lineage]}.vcf.gz"),
                                  vcfFileIndex = glue::glue("{G2G_Obj$host_path[cur_lineage]}.vcf.gz.tbi"),
                                  vcfField = "GT",
                                  chrom = as.character(j),
                                  minMAC = 4,
                                  IsDropMissingDosages = T,
                                  IsOutputAFinCaseCtrl = T,
                                  IsOutputNinCaseCtrl = T,
                                  minMAF = 0.05,
                                  LOCO = F,
                                  numLinesOutput = .Machine$integer.max,
                                  sampleFile= glue::glue("{G2G_Obj$host_path[cur_lineage]}.vcf.gz.sample"),
                                  GMMATmodelFile= glue::glue('{OUT_PATH_Lineage}{cur_pathogen_variant}.rda'),
                                  varianceRatioFile=glue::glue('{OUT_PATH_Lineage}{cur_pathogen_variant}.varianceRatio.txt'),
                                  SAIGEOutputFile= glue::glue('{OUT_PATH_Lineage}/SAIGE_SPA_TEST_chr_{j}_{cur_pathogen_variant}'))
             T}
              ,error = function(e){
                write(glue::glue("CHR {j} ERROR"),glue::glue('{OUT_PATH_Lineage}/SAIGE_SPA_TEST_{cur_pathogen_variant}.error'),append = T)
                write(as.character(e),glue::glue('{{OUT_PATH_Lineage}/SAIGE_SPA_TEST_{cur_pathogen_variant}.error'),append = T)
                return(F)},
              warning = function(w){
                write(glue::glue("CHR {j} WARNING"),glue::glue('{OUT_PATH_Lineage}/SAIGE_SPA_TEST_{cur_pathogen_variant}.warning'),append = T)
                write(as.character(w),glue::glue('{OUT_PATH_Lineage}/SAIGE_SPA_TEST_{cur_pathogen_variant}.warning'),append = T)
                return(F)})
            return(SPA_test_stat)
          },mc.preschedule = F,mc.cores = n_cores)
        }else if (tool == 'GMMAT-WALD' | tool == 'GMMAT-SCORE'){
          gmmat_pheno <- dplyr::left_join(data.frame(cur_aa_Matrix_mat),covars,by=c('IID'='IID','FID'='FID')) %>% dplyr::select(-FID) %>% dplyr::rename(ID = IID)
          grm_mat <- read.table(file = glue::glue("{OUT_PATH_Lineage}/tmp/output/GEMMA.sXX.txt") ,check.names = F)
          grm_mat <- as.matrix(grm_mat)
          fam_file <- data.table::fread(glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"),header = F)
          colnames(grm_mat) <- fam_file$V2
          rownames(grm_mat) <- fam_file$V2
          
          
          covar_names <- paste0(setdiff(colnames(covars),c('FID','IID')),collapse = '+')
          if(tool == 'GMMAT-SCORE'){
            null_model <- GMMAT::glmmkin(glue::glue("{cur_pathogen_variant} ~ {covar_names}"),
                                         data = gmmat_pheno,
                                         kins = grm_mat,
                                         id = "ID",
                                         family = binomial(link = "logit"))
            GMMAT::glmm.score(null_model,
                              infile = glue::glue('{G2G_Obj$host_path[cur_lineage]}'),
                              outfile = glue::glue('{OUT_PATH_Lineage}/GMMAT_score_{cur_pathogen_variant}.pval'))
            system(glue::glue('pigz --fast {OUT_PATH_Lineage}/GMMAT_score_{cur_pathogen_variant}.pval'))
            
          }else if(tool == 'GMMAT-WALD'){
            wald_pheno <- dplyr::select(gmmat_pheno,ID,AA = cur_pathogen_variant,PC1,PC2,PC3)
            wald_test_result <- GMMAT::glmm.wald(fixed = AA ~ PC1+PC2+PC3,
                                                 data = wald_pheno,
                                                 kins = grm_mat,
                                                 id = "ID",
                                                 family = binomial(link = "logit"),
                                                 infile = glue::glue("{G2G_Obj$host_path[cur_lineage]}"),
                                                 snps = test_snp,
                                                 verbose=T)
            data.table::fwrite(wald_test_result,glue::glue('{OUT_PATH_Lineage}//GMMAT_wald_{cur_pathogen_variant}.pval'),sep = ' ',col.names = T,row.names = F)
          }
        }
      }
    }
  }
}

RunG2GPerm <- function(G2G_Obj,SOFTWARE_DIR,OUT_DIR,Ref_Panel,tool = 'PLINK',n_PC = 5,n_pPC = 6,n_cores = 20,covars_to_incl = c(),lineage = c(),debug=F,chunks = c(),n_Perm = 1){
  OUT_PATH <- glue::glue("{OUT_DIR}/{Ref_Panel}")
  system(glue::glue("mkdir -p {OUT_PATH}"))
  
  PLINK = '/home/zmxu/Software/plink'
  bcftools = '/home/zmxu/Software/bcftools'
  OUT_PATH <- glue::glue("{OUT_PATH}/{tool}/")
  system(glue::glue("mkdir -p {OUT_PATH}"))
  if(length(covars_to_incl) > 0){
    OUT_PATH <- glue::glue("{OUT_PATH}/PC_{n_PC}_pPC_{n_pPC}_cov_{paste0(sapply(covars_to_incl,function(x) gsub(x=x,pattern='_',replacement='')),collapse = '-')}/")
  }else{
    OUT_PATH <- glue::glue("{OUT_PATH}/PC_{n_PC}_pPC_{n_pPC}/")
  }
  system(glue::glue("mkdir -p {OUT_PATH}"))
  
  if (tool == 'PLINK' | tool == 'PLINK-FIRTH' | tool == 'GMMAT-WALD' | tool == 'GMMAT-SCORE'){
    if(length(lineage)==0){
      lineages_to_run <- unique(names(G2G_Obj$aa_matrix_filt))
    }else{
      lineages_to_run <- lineage
    }
    for(i in 1:length(lineages_to_run)){
      cur_lineage <- lineages_to_run[i]
      OUT_PATH_Lineage <- glue::glue("{OUT_PATH}/LINEAGE_{cur_lineage}/")
      system(glue::glue("mkdir -p {OUT_PATH_Lineage}"))
      system(glue::glue("mkdir -p {OUT_PATH_Lineage}/tmp"))
      
      #Get IDD and FID, ensure they're in correct order
      fam_file <- data.table::fread(glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"))
      colnames(fam_file)[1:2] <- c('FID','IID')
      if(!all(fam_file$IID == G2G_Obj$both_IDs_to_keep[[cur_lineage]]$FAM_ID)){
        stop('Sample order incorrect')
      }
      
      #Write out AA matrix
      cur_aa_matrix_filt <- G2G_Obj$aa_matrix_filt[[cur_lineage]]
      cur_aa_Matrix_mat <- cbind(fam_file[,c(1,2)],as.matrix(cur_aa_matrix_filt))
      colnames(cur_aa_Matrix_mat) <- c('FID','IID',colnames(cur_aa_matrix_filt))
      
      if(debug){
        cur_aa_Matrix_mat <- cur_aa_Matrix_mat[,1:5]
        data.table::fwrite(cur_aa_Matrix_mat,col.names = T,row.names = F,sep = ' ',file = glue::glue("{OUT_PATH_Lineage}/tmp/AA_outcome.txt"),na = 'NA',quote = F)
      }else{
        data.table::fwrite(cur_aa_Matrix_mat,col.names = T,row.names = F,sep = ' ',file = glue::glue("{OUT_PATH_Lineage}/tmp/AA_outcome.txt"),na = 'NA',quote = F)
      }
      
      #Store PCs and covars
      host_PCs <- G2G_Obj$host_PCs[[cur_lineage]]
      host_PCs <- dplyr::select(host_PCs,'IID',paste0('PC',1:n_PC))
      
      if(n_pPC != 0){
        pPCs <- G2G_Obj$vir_pPCs[[cur_lineage]]
        pPCs <-dplyr::select(pPCs,IID,paste0('PC',1:n_pPC))
        colnames(pPCs) <- c('IID',paste0('pPC',1:n_pPC))
      }
      if(length(covars_to_incl) > 0){
        covars_num_discrete <- G2G_Obj$covars[[cur_lineage]]
        covars_num_discrete <- dplyr::select(covars_num_discrete,'IID',covars_to_incl)
        covars_num_discrete[is.na(covars_num_discrete)] <- 'NONE'
        
        if(n_pPC != 0){
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(pPCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(covars_num_discrete,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }else{
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(covars_num_discrete,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }
      }else{
        if(n_pPC != 0){
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(pPCs,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }else if(n_PC != 0){
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }else{
          covars <- fam_file %>% dplyr::select(FID,IID)
        }
      }
      
      #Change categorial binary covariates to numeric (0 and 1)
      if(ncol(covars) > 2){
        cat_covars <- which(sapply(3:ncol(covars),function(q) any(is.character(as.vector(t(covars[,..q]))))))
        if(length(cat_covars) > 0){
          for(col in (cat_covars+2)){
            covars[,col] <- as.integer(as.factor(t(covars[,..col]))) - 1
          }
        }
      }

      data.table::fwrite(covars,col.names = F,row.names = F,sep = ' ',file = glue::glue("{OUT_PATH_Lineage}/tmp/plink-covars.txt"),na = 'NA',quote = F)
      
      #Run association study for each AA variant in the current lineage
      AA_Matrix_No_ID <- cur_aa_Matrix_mat[,-c(1,2)]
      #Run association study for each AA variant in the current lineage
      if(length(chunks) == 0){
        chunks <- 1:ncol(AA_Matrix_No_ID)
      }
      if(max(chunks) > ncol(AA_Matrix_No_ID)){
        chunks <- chunks[1]:ncol(AA_Matrix_No_ID)
      }
      
      #Generate Permutations
      if(file.exists(glue::glue("{OUT_PATH_Lineage}/tmp/perm_list.rds"))){
        perm_list <- readRDS(glue::glue("{OUT_PATH_Lineage}/tmp/perm_list.rds"))
      }else{
        perm_list <- lapply(1:n_Perm,function(x) sample(1:nrow(fam_file),size = nrow(fam_file),replace = F))
        saveRDS(perm_list,file = glue::glue("{OUT_PATH_Lineage}/tmp/perm_list.rds"))
      }
      #Preprocess for GMATT
      if(tool == 'GMMAT-WALD' | tool == 'GMMAT-SCORE'){
        if(!file.exists(glue::glue("{OUT_PATH_Lineage}/tmp/output/GEMMA.sXX.txt"))){
          #Run GEMMA to make GRM (https://www.xzlab.org/software/GEMMAmanual.pdf)
          cur_dir <- getwd()
          setwd(glue::glue("{OUT_PATH_Lineage}/tmp/"))
          fam_file <- data.table::fread(glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"),header = F)
          fam_file$V6 <- 0 #GEMMA dislikes -9 label (but doesn't use pheno info)
          data.table::fwrite(fam_file,glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"),col.names = F,sep = ' ',append = F)
          
          system(glue::glue("{SOFTWARE_DIR}/gemma-0.98.1-linux-static -bfile {G2G_Obj$host_path[cur_lineage]} -gk 2 -o /GEMMA"))
          setwd(cur_dir)
        }
        gmmat_pheno <- dplyr::left_join(data.frame(cur_aa_Matrix_mat),covars,by=c('IID'='IID','FID'='FID')) %>% dplyr::select(-FID) %>% dplyr::rename(ID = IID)
        grm_mat <- read.table(file = glue::glue("{OUT_PATH_Lineage}/tmp/output/GEMMA.sXX.txt") ,check.names = F)
        grm_mat <- as.matrix(grm_mat)
        fam_file <- data.table::fread(glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"),header = F)
        colnames(grm_mat) <- fam_file$V2
        rownames(grm_mat) <- fam_file$V2
      }
      
      for(q in 1:length(perm_list)){
        if(tool == 'PLINK' & !file.exists(glue::glue("{OUT_PATH_Lineage}/tmp/perm{q}.fam"))){
          cur_perm_fam <- fam_file[perm_list[[q]],]
          data.table::fwrite(cur_perm_fam,sep = ' ',col.names = F,row.names = F,quote = F,file = glue::glue("{OUT_PATH_Lineage}/tmp/perm{q}.fam"))
          system(glue::glue("mkdir -p {OUT_PATH_Lineage}/perm{q}/"))
        }
        for(k in chunks){
          cur_pathogen_variant <- colnames(AA_Matrix_No_ID)[k]
          if(tool == 'PLINK'){
            if(ncol(covars) > 2){
              tryCatch(system(
                glue::glue(
                  "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --fam {OUT_PATH_Lineage}/tmp/perm{q}.fam --no-sex --glm cc-residualize hide-covar --covar-variance-standardize --1 --pheno {OUT_PATH_Lineage}/tmp/AA_outcome.txt --pheno-col-nums {k+2} --covar {OUT_PATH_Lineage}/tmp/plink-covars.txt --out {OUT_PATH_Lineage}/perm{q}/{cur_pathogen_variant}"
                )
              ))
            }else{
              tryCatch(system(
                glue::glue(
                  "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --fam {OUT_PATH_Lineage}/tmp/perm{q}.fam --no-sex --glm cc-residualize hide-covar allow-no-covars --1 --pheno {OUT_PATH_Lineage}/tmp/AA_outcome.txt --pheno-col-nums {k+2} --out {OUT_PATH_Lineage}/perm{q}/{cur_pathogen_variant}"
                )
              ))
            }
            system(glue::glue("awk \'{{print $13\" \"$14}}\' {OUT_PATH_Lineage}/perm{q}/{cur_pathogen_variant}*glm.logistic.hybrid > {OUT_PATH_Lineage}/perm{q}/{cur_pathogen_variant}.pval"))
            system(glue::glue('pigz --best {OUT_PATH_Lineage}/perm{q}/{cur_pathogen_variant}.pval'))
          }else if(tool == 'GMMAT-SCORE'){
            system(glue::glue("mkdir -p {OUT_PATH_Lineage}/perm{q}/"))
            
            cur_perm_pheno <- data.frame(ID = gmmat_pheno[,'ID'],gmmat_pheno[perm_list[[q]],!colnames(gmmat_pheno)=='ID'])
            covar_names <- paste0(setdiff(colnames(covars),c('FID','IID')),collapse = '+')
            null_model <- GMMAT::glmmkin(glue::glue("{cur_pathogen_variant} ~ {covar_names}"),
                                         data = cur_perm_pheno,
                                         kins = grm_mat,
                                         id = "ID",
                                         family = binomial(link = "logit"))
            GMMAT::glmm.score(null_model,
                              infile = glue::glue('{G2G_Obj$host_path[cur_lineage]}'),
                              outfile = glue::glue('{OUT_PATH_Lineage}/perm{q}/{cur_pathogen_variant}.raw'))
            system(glue::glue("awk  '{{print $11}}' {OUT_PATH_Lineage}/perm{q}/{cur_pathogen_variant}.raw > {OUT_PATH_Lineage}/perm{q}/{cur_pathogen_variant}.pval"))
            
            system(glue::glue('pigz --fast {OUT_PATH_Lineage}/perm{q}/{cur_pathogen_variant}.pval'))
            system(glue::glue("rm {OUT_PATH_Lineage}/perm{q}/{cur_pathogen_variant}.raw"))
          }
        }
      }
    }
  }
}

RunInteraction <- function(G2G_Obj,SOFTWARE_DIR,OUT_DIR,Ref_Panel,tool = 'PLINK',n_PC = 3,n_pPC = 3,n_cores = 20,covars_to_incl = c(),lineage = c(),by_lineage = F,chunks = c(),debug=F){
  OUT_PATH <- glue::glue("{OUT_DIR}/{Ref_Panel}")
  system(glue::glue("mkdir -p {OUT_PATH}"))
  
  PLINK = '/home/zmxu/Software/plink'
  GCTA = '/home/zmxu/Software/gcta'
  
  OUT_PATH <- glue::glue("{OUT_PATH}/{tool}/")
  system(glue::glue("mkdir -p {OUT_PATH}"))
  if(length(covars_to_incl) > 0){
    OUT_PATH <- glue::glue("{OUT_PATH}/PC_{n_PC}_pPC_{n_pPC}_cov_{paste0(sapply(covars_to_incl,function(x) gsub(x=x,pattern='_',replacement='')),collapse = '-')}/")
  }else{
    OUT_PATH <- glue::glue("{OUT_PATH}/PC_{n_PC}_pPC_{n_pPC}/")
  }
  system(glue::glue("mkdir -p {OUT_PATH}"))
  
  if (tool == 'PLINK'){
    if(length(lineage)==0){
      lineages_to_run <- unique(names(G2G_Obj$both_IDs_to_keep))
    }else{
      lineages_to_run <- lineage
    }
    for(i in 1:length(lineages_to_run)){
      cur_lineage <- lineages_to_run[i]
      OUT_PATH_Lineage <- glue::glue("{OUT_PATH}/LINEAGE_{cur_lineage}/")
      system(glue::glue("mkdir -p {OUT_PATH_Lineage}"))
      system(glue::glue("mkdir -p {OUT_PATH_Lineage}/tmp"))
      
      #Get IDD and FID, ensure they're in correct order
      fam_file <- data.table::fread(glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"))
      colnames(fam_file)[1:2] <- c('FID','IID')
      if(!all(fam_file$IID == G2G_Obj$both_IDs_to_keep[[cur_lineage]]$FAM_ID)){
        stop('Sample order incorrect')
      }
      
      #Store AA matrix
      if(cur_lineage == 'ALL' & !by_lineage){
        cur_aa_matrix_filt <- G2G_Obj$aa_matrix_raw
      }else if(cur_lineage == 'ALL' & by_lineage){
        #SNP x LINEAGE interaction
        cur_aa_matrix_filt <- as.matrix(one_hot(as.data.table(data.frame(LINEAGE = G2G_Obj$vir_pPCs$ALL$LINEAGE))))
      }else{
        cur_aa_matrix_filt <- G2G_Obj$aa_matrix_filt[[cur_lineage]]
      }
      
      #Store PCs and covars
      host_PCs <- G2G_Obj$host_PCs[[cur_lineage]]
      host_PCs <- dplyr::select(host_PCs,'IID',paste0('PC',1:n_PC))
      
      if(n_pPC != 0){
        pPCs <- G2G_Obj$vir_pPCs[[cur_lineage]]
        pPCs <-dplyr::select(pPCs,IID,paste0('PC',1:n_pPC))
        colnames(pPCs) <- c('IID',paste0('pPC',1:n_pPC))
      }
      if(length(covars_to_incl) > 0){
        covars_num_discrete <- G2G_Obj$covars[[cur_lineage]]
        covars_num_discrete <- dplyr::select(covars_num_discrete,'IID',covars_to_incl)
        if(n_pPC != 0){
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(pPCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(covars_num_discrete,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }else{
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(covars_num_discrete,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }
      }else{
        if(n_pPC != 0){
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
            dplyr::left_join(pPCs,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }else{
          covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
        }
      }
      #Change categorial binary covariates to numeric (0 and 1)
      cat_covars <- which(sapply(3:ncol(covars),function(q) any(is.character(as.vector(t(covars[,..q]))))))
      if(length(cat_covars) > 0){
        for(col in (cat_covars+2)){
          covars[,col] <- as.integer(as.factor(t(covars[,..col]))) - 1
        }
      }
      
      #Write out Phenotype
      tb_score <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),G2G_Obj$tb_score[[cur_lineage]] %>% dplyr::select(IID,tb_score),by=c('IID'='IID')) %>% dplyr::relocate(FID,IID,tb_score)
      data.table::fwrite(tb_score,col.names = T,quote = F,row.names = F,sep = ' ',file = glue::glue("{OUT_PATH_Lineage}/tmp/tb_score.txt"))
      
      #Run association study for each AA variant in the current lineage
      if(length(chunks) == 0){
        chunks <- 1:ncol(cur_aa_matrix_filt)
      }
      if(max(chunks) > ncol(cur_aa_matrix_filt)){
        chunks <- chunks[1]:ncol(cur_aa_matrix_filt)
      }
      #Run GWAS with lineage as covariate
      if(by_lineage){
        cur_covar <- cbind(covars[,1:2],data.frame(LINEAGE = G2G_Obj$vir_pPCs$ALL$LINEAGE),covars[,-c(1,2)])
        data.table::fwrite(cur_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_PATH_Lineage}/tmp/lineage_covar.txt"))

        system(
          glue::glue(
            "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --covar-variance-standardize --linear no-x-sex --pheno {OUT_PATH_Lineage}/tmp/tb_score.txt --covar {OUT_PATH_Lineage}/tmp/lineage_covar.txt --out {OUT_PATH_Lineage}Lineage_GWAS"
          )
        )

        discrete_covar <- dplyr::select(cur_covar,c('FID','IID',intersect(colnames(cur_covar),c('Patient_Sex','HIV_Status','LINEAGE'))))
        num_covar <- dplyr::select(cur_covar,c('FID','IID',setdiff(colnames(cur_covar),colnames(discrete_covar))))
        data.table::fwrite(discrete_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_PATH_Lineage}/tmp/lineage_discrete_covar.txt"))
        data.table::fwrite(num_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_PATH_Lineage}/tmp/lineage_num_covar.txt"))
        system(
          glue::glue(
            "{GCTA} --bfile {G2G_Obj$host_path[cur_lineage]} --make-grm --out {G2G_Obj$host_path[cur_lineage]} --thread-num {n_cores}"
          )
        )

        system(
          glue::glue(
            "{GCTA} --mlma --bfile {G2G_Obj$host_path[cur_lineage]} --grm {G2G_Obj$host_path[cur_lineage]} --autosome --pheno {OUT_PATH_Lineage}/tmp/tb_score.txt --covar {OUT_PATH_Lineage}/tmp/lineage_discrete_covar.txt --qcovar {OUT_PATH_Lineage}/tmp/lineage_num_covar.txt --thread-num {n_cores} --out {OUT_PATH_Lineage}Lineage_GWAS"
          )
        )

      }
      
      
      for(k in chunks){
        tryCatch({
          if(by_lineage){
            cur_pathogen_variant <- colnames(cur_aa_matrix_filt)[k]
            cur_covar <- cbind(covars[,1:2],cur_aa_matrix_filt[,k,drop=FALSE],cur_aa_matrix_filt[,-k][,-1],covars[,-c(1,2)])
            data.table::fwrite(cur_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_PATH_Lineage}/tmp/{cur_pathogen_variant}.txt"))
            
            system(
              glue::glue(
                "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --covar-variance-standardize --linear no-x-sex interaction --parameters 1-{ncol(cur_covar)} --pheno {OUT_PATH_Lineage}/tmp/tb_score.txt --covar {OUT_PATH_Lineage}/tmp/{cur_pathogen_variant}.txt --out {OUT_PATH_Lineage}{cur_pathogen_variant}"
              )
            )
          }else{
            cur_pathogen_variant <- colnames(cur_aa_matrix_filt)[k]
            cur_covar <- cbind(covars[,1:2],cur_aa_matrix_filt[,k,drop=FALSE],covars[,-c(1,2)])
            data.table::fwrite(cur_covar,col.names = T,quote = F,row.names = F,sep = ' ',na = 'NA',file = glue::glue("{OUT_PATH_Lineage}/tmp/{cur_pathogen_variant}.txt"))
            
            system(
              glue::glue(
                "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --covar-variance-standardize --linear no-x-sex interaction --parameters 1-{ncol(cur_covar)} --pheno {OUT_PATH_Lineage}/tmp/tb_score.txt --covar {OUT_PATH_Lineage}/tmp/{cur_pathogen_variant}.txt --out {OUT_PATH_Lineage}{cur_pathogen_variant}"
              )
            )
            system(glue::glue("grep 'ADD' {OUT_PATH_Lineage}{cur_pathogen_variant}.tb_score.glm.linear > {OUT_PATH_Lineage}{cur_pathogen_variant}.tb_score.glm.linear.add"))
            system(glue::glue("{SOFTWARE_DIR}pigz --fast {OUT_PATH_Lineage}{cur_pathogen_variant}.tb_score.glm.linear.add"))
            system(glue::glue("rm {OUT_PATH_Lineage}{cur_pathogen_variant}.tb_score.glm.linear"))
          }
        }
        )
      }
      # results <- GetResults(OUT_PATH_Lineage,suffix = 'glm.linear.add.gz',p_thresh = 5e-8,n_cores = n_cores,is_interaction = T,is_ordinal = F)
      # saveRDS(results,glue::glue("{OUT_PATH_Lineage}results.rds"))
      
    }
  }else if (tool == 'ORDINAL' | tool == 'ORDINAL_GxE'){
    if(length(lineage)==0){
      lineages_to_run <- unique(names(G2G_Obj$both_IDs_to_keep))
    }else{
      lineages_to_run <- lineage
    }
    for(i in 1:length(lineages_to_run)){
      cur_lineage <- lineages_to_run[i]
      OUT_PATH_Lineage <- glue::glue("{OUT_PATH}/LINEAGE_{cur_lineage}/")
      system(glue::glue("mkdir -p {OUT_PATH_Lineage}"))
      system(glue::glue("mkdir -p {OUT_PATH_Lineage}/tmp"))
      
      cur_wd <- getwd()
      setwd(OUT_PATH_Lineage)
      
      
      #Get IDD and FID, ensure they're in correct order
      fam_file <- data.table::fread(glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"))
      colnames(fam_file)[1:2] <- c('FID','IID')
      if(!all(fam_file$IID == G2G_Obj$both_IDs_to_keep[[cur_lineage]]$FAM_ID)){
        stop('Sample order incorrect')
      }
      
      #Store AA matrix
      if(cur_lineage == 'ALL' & !by_lineage){
        cur_aa_matrix_filt <- G2G_Obj$aa_matrix_raw
      }else if(cur_lineage == 'ALL' & by_lineage){
        #SNP x LINEAGE interaction
        cur_aa_matrix_filt <- as.matrix(one_hot(as.data.table(data.frame(LINEAGE = G2G_Obj$vir_pPCs$ALL$LINEAGE))))
      }else{
        cur_aa_matrix_filt <- G2G_Obj$aa_matrix_filt[[cur_lineage]]
      }
      
      #Store PCs and covars
      host_PCs <- G2G_Obj$host_PCs[[cur_lineage]]
      host_PCs <- dplyr::select(host_PCs,paste0('PC',1:n_PC))
      
      if(n_pPC != 0){
        pPCs <- G2G_Obj$vir_pPCs[[cur_lineage]]
        pPCs <-dplyr::select(pPCs,paste0('PC',1:n_pPC))
        colnames(pPCs) <- paste0('pPC',1:n_pPC)
      }
      if(length(covars_to_incl) > 0){
        covars_num_discrete <- G2G_Obj$covars[[cur_lineage]]
        covars_num_discrete <- dplyr::select(covars_num_discrete,covars_to_incl)
        if(n_pPC != 0){
          covars <- cbind(fam_file[,c(1,2)],host_PCs,pPCs,covars_num_discrete)
        }else{
          covars <- cbind(fam_file[,c(1,2)],host_PCs,covars_num_discrete)
        }
      }else{
        if(n_pPC != 0){
          covars <- cbind(fam_file[,c(1,2)],host_PCs,pPCs)
        }else{
          covars <- cbind(fam_file[,c(1,2)],host_PCs)
        }
      }
      
      #Store Phenotype
      pheno <- G2G_Obj$tb_score[[cur_lineage]]
      
      outcome_name <- setdiff(colnames(pheno),c('FID','IID'))
      covar_names <- paste(setdiff(colnames(covars),c('FID','IID')),collapse = ' + ')
      
      #Split genotype by chromosome
      lapply(1:22,function(q) system(glue::glue("{PLINK}2 --bfile {G2G_Obj$host_path[cur_lineage]} --set-all-var-ids @:#[b37]\\$r_\\$a --chr {q} --make-bed --out {OUT_PATH_Lineage}/tmp/chr{q}")))
      
      #Run association study for each AA variant in the current lineage
      if(length(chunks) == 0){
        chunks <- 1:ncol(cur_aa_matrix_filt)
      }
      if(max(chunks) > ncol(cur_aa_matrix_filt)){
        chunks <- chunks[1]:ncol(cur_aa_matrix_filt)
      }
      for(k in chunks){
        print(k)
        cur_pathogen_variant <- colnames(cur_aa_matrix_filt)[k]
        cur_pathogen_dosage <- data.frame(AA = cur_aa_matrix_filt[,k,drop=T])
        cur_pheno <- cbind(data.frame(famid=covars[,1],perid=covars[,2],faid=0,moid=0),cur_pathogen_dosage,covars[,-c(1,2)],pheno[,-c(1,2)])
        write.csv(cur_pheno,quote = F,row.names = F,file = glue::glue("{OUT_PATH_Lineage}/tmp/{cur_pathogen_variant}.txt"))
        
        tryCatch({
          if(tool == 'ORDINAL'){
            mclapply(1:22,function(q){
              system(glue::glue("nice -n 19 julia -e 'using OrdinalGWAS; const datadir = \"{OUT_PATH_Lineage}/tmp/\"; ordinalgwas(@formula({outcome_name} ~ AA + {covar_names}), datadir * \"{cur_pathogen_variant}.txt\", datadir * \"chr{q}\",pvalfile=\"chr{q}\",testformula=@formula({outcome_name} ~ snp + snp & AA))'"))
            },mc.cores = n_cores,mc.silent = T,mc.preschedule = F)
            
          }else if(tool == 'ORDINAL_GxE'){
            mclapply(1:22,function(q){
              system(glue::glue("nice -n 19 julia -e 'using OrdinalGWAS; const datadir = \"{OUT_PATH_Lineage}/tmp/\"; ordinalgwas(@formula({outcome_name} ~ AA + {covar_names}), datadir * \"{cur_pathogen_variant}.txt\", datadir * \"chr{q}\",pvalfile=\"chr{q}\",analysistype = \"gxe\",e = :AA,test = :score)'"))
            },mc.cores = n_cores,mc.silent = T,mc.preschedule = F)
          }
          system(glue::glue("head -n 1 {OUT_PATH_Lineage}chr1 > {OUT_PATH_Lineage}{cur_pathogen_variant}.pval"))
          trash <- sapply(1:22,function(q) system(glue::glue("cat {OUT_PATH_Lineage}chr{q} | grep -v 'chr' >> {OUT_PATH_Lineage}{cur_pathogen_variant}.pval")))
          system(glue::glue("mv {OUT_PATH_Lineage}ordinalgwas.null.txt {OUT_PATH_Lineage}{cur_pathogen_variant}.null"))
          system(glue::glue("rm {OUT_PATH_Lineage}chr*"))
          system(glue::glue("pigz --fast {OUT_PATH_Lineage}{cur_pathogen_variant}.pval"))
        })
      }
      results <- GetResults(OUT_PATH_Lineage,suffix = '.pval.gz',p_thresh = 5e-8,n_cores = n_cores,is_interaction = T,is_ordinal = T)
      saveRDS(results,glue::glue("{OUT_PATH_Lineage}results.rds"))
      
    }
  }
}

RunG2GBurden <- function(G2G_Obj,SOFTWARE_DIR,OUT_DIR,Ref_Panel,tool = 'PLINK',n_PC = 3,n_pPC = 0,n_cores = 20,covars_to_incl = c(),gene = c(),X_Chr=T,MAC_Thresh = 10,binarize = T,incl_syn = F){
  OUT_PATH <- glue::glue("{OUT_DIR}/{Ref_Panel}")
  system(glue::glue("mkdir -p {OUT_PATH}"))
  
  PLINK = '/home/zmxu/Software/plink'
  
  OUT_PATH <- glue::glue("{OUT_PATH}/{tool}/")
  system(glue::glue("mkdir -p {OUT_PATH}"))
  if(length(covars_to_incl) > 0 & incl_syn){
    OUT_PATH <- glue::glue("{OUT_PATH}/PC_{n_PC}_pPC_{n_pPC}_cov_{paste0(sapply(covars_to_incl,function(x) gsub(x=x,pattern='_',replacement='')),collapse = '-')}_Syn/")
  }else if(length(covars_to_incl) > 0){
    OUT_PATH <- glue::glue("{OUT_PATH}/PC_{n_PC}_pPC_{n_pPC}_cov_{paste0(sapply(covars_to_incl,function(x) gsub(x=x,pattern='_',replacement='')),collapse = '-')}/")
  }else if(incl_syn){
    OUT_PATH <- glue::glue("{OUT_PATH}/PC_{n_PC}_pPC_{n_pPC}_Syn/")
  }else{
    OUT_PATH <- glue::glue("{OUT_PATH}/PC_{n_PC}_pPC_{n_pPC}/")
  }
  system(glue::glue("mkdir -p {OUT_PATH}"))
  
  if (tool == 'PLINK'){
    cur_lineage <- 'ALL'
    OUT_PATH_Lineage <- glue::glue("{OUT_PATH}/LINEAGE_{cur_lineage}/")
    system(glue::glue("mkdir -p {OUT_PATH_Lineage}"))
    system(glue::glue("mkdir -p {OUT_PATH_Lineage}/tmp"))
    
    #Get IDD and FID, ensure they're in correct order
    fam_file <- data.table::fread(glue::glue("{G2G_Obj$host_path[cur_lineage]}.fam"))
    colnames(fam_file)[1:2] <- c('FID','IID')
    if(!all(fam_file$IID == G2G_Obj$both_IDs_to_keep[[cur_lineage]]$FAM_ID)){
      stop('Sample order incorrect')
    }
    
    #Write out Burden Matrix
    cur_burden_matrix <- G2G_Obj$gene_burden_non_syn
    MAC <- apply(cur_burden_matrix,2,function(x) min(sum(x!=0),sum(x==0)))
    cur_burden_matrix_filt <- cur_burden_matrix[,MAC > MAC_Thresh]
    #Binarize matrix if specified (dominant, case-control model)
    if(binarize){
      cur_burden_matrix_filt <- apply(cur_burden_matrix_filt,2,function(x) ifelse(x!=0,1,0))
    }
    
    cur_burden_Matrix_mat <- cbind(fam_file[,c(1,2)],as.matrix(cur_burden_matrix_filt))
    colnames(cur_burden_Matrix_mat) <- c('FID','IID',colnames(cur_burden_matrix_filt))
    
    #Write out burden matrix
    data.table::fwrite(cur_burden_Matrix_mat,col.names = T,row.names = F,sep = ' ',file = glue::glue("{OUT_PATH_Lineage}/tmp/Burden_outcome.txt"),na = 'NA',quote = F)
    
    #Store PCs and covars
    host_PCs <- G2G_Obj$host_PCs[[cur_lineage]]
    host_PCs <- dplyr::select(host_PCs,'IID',paste0('PC',1:n_PC))
    
    if(n_pPC != 0){
      pPCs <- G2G_Obj$vir_pPCs[[cur_lineage]]
      pPCs <- dplyr::select(pPCs,PATIENT_ID,paste0('PC',1:n_pPC))
      colnames(pPCs) <- c('PATIENT_ID',paste0('pPC',1:n_pPC))
    }
    if(length(covars_to_incl) > 0){
      covars_num_discrete <- G2G_Obj$covars[[cur_lineage]]
      covars_num_discrete <- dplyr::select(covars_num_discrete,'IID',covars_to_incl)
      covars_num_discrete[is.na(covars_num_discrete)] <- 'NONE'
      
      if(n_pPC != 0){
        covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
          dplyr::left_join(pPCs %>% dplyr::left_join(G2G_Obj$both_IDs_to_keep[[cur_lineage]] %>% dplyr::select(PATIENT_ID,IID=FAM_ID),by=c('PATIENT_ID' = 'PATIENT_ID')) %>% dplyr::select(-PATIENT_ID),by=c('IID'='IID')) %>% 
          dplyr::left_join(covars_num_discrete,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
      }else{
        covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
          dplyr::left_join(covars_num_discrete,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
      }
    }else{
      if(n_pPC != 0){
        covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% 
          dplyr::left_join(pPCs %>% dplyr::left_join(G2G_Obj$both_IDs_to_keep[[cur_lineage]] %>% dplyr::select(PATIENT_ID,IID=FAM_ID),by=c('PATIENT_ID' = 'PATIENT_ID')) %>% dplyr::select(-PATIENT_ID),by=c('IID'='IID')) %>%
          dplyr::relocate(FID,IID)
      }else{
        covars <- dplyr::left_join(fam_file %>% dplyr::select(FID,IID),host_PCs,by=c('IID'='IID')) %>% dplyr::relocate(FID,IID)
      }
    }
    
    #Change categorial binary covariates to numeric (0 and 1)
    cat_covars <- which(sapply(3:ncol(covars),function(q) any(is.character(as.vector(t(covars[,..q]))))))
    if(length(cat_covars) > 0){
      for(col in (cat_covars+2)){
        covars[,col] <- as.integer(as.factor(t(covars[,..col]))) - 1
      }
    }
    data.table::fwrite(covars,col.names = F,row.names = F,sep = ' ',file = glue::glue("{OUT_PATH_Lineage}/tmp/plink-covars.txt"),na = 'NA',quote = F)
    covar_path <- glue::glue("{OUT_PATH_Lineage}/tmp/plink-covars.txt")
    
    #Run association study for each AA variant in the current lineage
    Burden_Matrix_No_ID <- cur_burden_Matrix_mat[,-c(1,2)]
    #Run association study for each AA variant in the current lineage
    if(length(gene) == 0 & !incl_syn){
      chunks <- 1:ncol(Burden_Matrix_No_ID)
    }else if(length(gene) == 0 & incl_syn){
      chunks <- match(intersect(colnames(Burden_Matrix_No_ID),colnames(G2G_Obj$gene_burden_syn)),colnames(Burden_Matrix_No_ID))
    }else{
      chunks <- match(gene,colnames(Burden_Matrix_No_ID))
    }
    genes <- colnames(Burden_Matrix_No_ID)
    for(k in 1:length(chunks)){
      cur_pathogen_gene <- genes[chunks[k]]
      var_std <- ' --covar-variance-standardize'
      if(incl_syn){
        cur_covars <- cbind(covars,data.frame(Syn_Count = G2G_Obj$gene_burden_syn[,match(cur_pathogen_gene,colnames(G2G_Obj$gene_burden_syn))]))
        data.table::fwrite(cur_covars,col.names = F,row.names = F,sep = ' ',file = glue::glue("{OUT_PATH_Lineage}/tmp/{cur_pathogen_gene}-covars.txt"),na = 'NA',quote = F)
        covar_path <- glue::glue("{OUT_PATH_Lineage}/tmp/{cur_pathogen_gene}-covars.txt")
        var_std <- ''
      }

      if(X_Chr){
        tryCatch(system(
          glue::glue(
            "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --glm cc-residualize hide-covar{var_std} --1 --pheno {OUT_PATH_Lineage}/tmp/Burden_outcome.txt --pheno-col-nums {chunks[k]+2} --covar {covar_path} --out {OUT_PATH_Lineage}{cur_pathogen_gene}"
          )
        ))

      }else{
        tryCatch(system(
          glue::glue(
            "{PLINK}2 --threads {n_cores} --bfile {G2G_Obj$host_path[cur_lineage]} --no-sex --glm cc-residualize hide-covar {var_std} --1 --pheno {OUT_PATH_Lineage}/tmp/Burden_outcome.txt --pheno-col-nums {chunks[k]+2} --covar {covar_path} --out {OUT_PATH_Lineage}{cur_pathogen_gene}"
          )
        ))

      }
      system(glue::glue('pigz --fast {OUT_PATH_Lineage}{cur_pathogen_gene}*.hybrid'))
    }
  }
}

ViewResults <- function(result_path,p_thresh,af_df=NULL,bonf=T,MAF_thresh = 0.05,is_interaction=F){
  results <- readRDS(result_path)
  if(is_interaction){
    results <- lapply(results,function(x) {
      if(ncol(x) == 0){
        return(data.frame('#CHROM'=character(),'POS'=numeric(),'ID'=character(),'REF'=character(),'ALT'=character(),'A1'=character(),'TEST'=character(),
                          'OBS_CT'=numeric(),'BETA'=numeric(),'SE'=numeric(),'T_STAT'=numeric(),'P'=numeric(),'ERRCODE'=character()))
      }else{
        colnames(x) <- c('#CHROM','POS','ID','REF','ALT','A1','TEST','OBS_CT','BETA','SE','T_STAT','P','ERRCODE')
        return(x)
      }
    })
  }
  if(!bonf){
    results_cleaned <- lapply(results,function(x) dplyr::filter(x,ERRCODE == '.',P < p_thresh))
  }else{
    results_cleaned <- lapply(results,function(x) dplyr::filter(x,ERRCODE == '.',P < p_thresh / length(results)))
  }
  results_cleaned_filt <- results_cleaned[sapply(results_cleaned,nrow) > 0]
  if(!is.null(af_df)){
    af_pass_thresh <- dplyr::filter(af_df,MAF > MAF_thresh)
    results_cleaned_jned <- mclapply(results_cleaned_filt,function(x) {
      dplyr::filter(x,ID %in% af_pass_thresh$SNP)
    },mc.cores = 10)
    results_cleaned_jned <- results_cleaned_jned[sapply(results_cleaned_jned,nrow) > 0]
    return(results_cleaned_jned)
  }else{
    return(results_cleaned_filt)
  }
}

SummaryStatsTbl <- function(G2G_Obj){
  
  metadata <- dplyr::left_join(G2G_Obj$both_IDs_to_keep[["ALL"]],G2G_Obj$tb_score$ALL %>% dplyr::select(IID,tb_score),by=c('FAM_ID' = 'IID'))
  metadata <- dplyr::left_join(metadata,by=c('FAM_ID' = 'IID'),G2G_Obj$covars$ALL %>% dplyr::select(-FID)) %>% 
    dplyr::left_join(G2G_Obj$host_PCs$ALL,by=c('FAM_ID' = 'IID')) %>%
    dplyr::left_join(G2G_Obj$vir_pPCs$ALL %>% dplyr::select(G_NUMBER,pPC1 = PC1,pPC2 = PC2, pPC3 = PC3),by=c('G_NUMBER' = 'G_NUMBER'))
  
  
  metadata$Patient_Sex[metadata$Patient_Sex == 1] <- 'female'
  metadata$Patient_Sex[metadata$Patient_Sex == 0] <- 'male'
  metadata$TB_RF_Smoking[metadata$TB_RF_Smoking == 1] <- 'yes'
  metadata$TB_RF_Smoking[metadata$TB_RF_Smoking == 0] <- 'no'
  
  lm("tb_score ~ Patient_Sex + Age + HIV_Status +LINEAGE + TB_RF_Smoking + PC1 + PC2 + PC3",data = metadata) %>% tbl_regression(tidy_fun = my_tidy) %>% bold_p()
  
  tbl_summary <- dplyr::select(metadata,LINEAGE,Age,PC1,PC2,PC3,TB_RF_Smoking,BMI,Patient_Household_Size)
  tbl_summary$Female <- ifelse(metadata$Patient_Sex == 'female',1,0)
  tbl_summary$HIV_Infected <- ifelse(metadata$HIV_Status == 'infected',1,0)
  tbl_summary(tbl_summary,by=LINEAGE,missing='no') %>% add_p() %>% add_overall() %>% bold_p()
  
} 


DATA_DIR <- '/home/zmxu/G2G_TB/data/'
SOFTWARE_DIR <- '/home/zmxu/Software/'

# G2G_Obj_AFGR <- SetUpG2G(DATA_DIR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/scratch/',Ref_Panel = 'AFGR',Host_MAF = 0.05,Pathogen_MAC_pPCA = 5,Pathogen_MAC_AA = 10)
# saveRDS(G2G_Obj_AFGR,'~/G2G_TB/scratch/AFGR/G2G_Obj_AFGR.rds')

# G2G_Obj_H3Africa <- SetUpG2G(DATA_DIR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/scratch/',Ref_Panel = 'H3Africa',Host_MAF = 0.05,Pathogen_MAC_pPCA = 5,Pathogen_MAC_AA = 10)
# saveRDS(G2G_Obj_H3Africa,'~/G2G_TB/scratch/H3Africa/G2G_Obj_H3Africa.rds')
# 
# G2G_Obj_Tanz <- SetUpG2G(DATA_DIR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/scratch/',Ref_Panel = 'Tanz',Host_MAF = 0.05,Pathogen_MAC_pPCA = 5,Pathogen_MAC_AA = 10)
# saveRDS(G2G_Obj_Tanz,'~/G2G_TB/scratch/Tanz/G2G_Obj_Tanz.rds')
# 
# G2G_Obj_TopMed <- SetUpG2G(DATA_DIR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/scratch/',Ref_Panel = 'TopMed',Host_MAF = 0.05,Pathogen_MAC_pPCA = 5,Pathogen_MAC_AA = 10)
# saveRDS(G2G_Obj_TopMed,'~/G2G_TB/scratch/TopMed/G2G_Obj_TopMed.rds')

# G2G_Obj_Tanz_AFGR <- SetUpG2G(DATA_DIR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/scratch/',Ref_Panel = 'Tanz_AFGR',Host_MAF = 0.05,Pathogen_MAC_pPCA = 5,Pathogen_MAC_AA = 10)
# saveRDS(G2G_Obj_Tanz_AFGR,'~/G2G_TB/scratch/Tanz_AFGR/G2G_Obj_Tanz_AFGR.rds')

# G2G_Obj_AFGR_Tanz <- SetUpG2G(DATA_DIR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/scratch/',Ref_Panel = 'AFGR_Tanz',Host_MAF = 0.05,Pathogen_MAC_pPCA = 5,Pathogen_MAC_AA = 10)
# saveRDS(G2G_Obj_AFGR_Tanz,'~/G2G_TB/scratch/AFGR_Tanz/G2G_Obj_AFGR_Tanz.rds')

# G2G_Obj_AFGR_Tanz_ChrX <- SetUpG2G(DATA_DIR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/scratch/',Ref_Panel = 'AFGR_Tanz_ChrX',Host_MAF = 0.05,Pathogen_MAC_pPCA = 5,Pathogen_MAC_AA = 10,VCF_Path = '~/TB_GWAS/data/Genotyping/TB_DAR_AFGR_Tanz_Imputed_ChrX/AFGR_Tanz.imputed.vcf.gz')
# saveRDS(G2G_Obj_AFGR_Tanz_ChrX,'~/G2G_TB/scratch/AFGR_Tanz_ChrX/G2G_Obj_AFGR_Tanz_ChrX.rds')

# G2G_Obj_AFGR_Tanz_ChrX_SIFT <- SetUpG2G(DATA_DIR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/scratch/',Ref_Panel = 'AFGR_Tanz_ChrX',Host_MAF = 0.05,Pathogen_MAC_pPCA = 5,Pathogen_MAC_AA = 10,VCF_Path = '~/TB_GWAS/data/Genotyping/TB_DAR_AFGR_Tanz_Imputed_ChrX/AFGR_Tanz.imputed.vcf.gz',SIFT = T)
# saveRDS(G2G_Obj_AFGR_Tanz_ChrX_SIFT,'~/G2G_TB/scratch/AFGR_Tanz_ChrX/G2G_Obj_AFGR_Tanz_ChrX_SIFT.rds')

# G2G_Obj_AFGR <- readRDS('~/G2G_TB/scratch/AFGR/G2G_Obj_AFGR.rds')
# G2G_Obj_H3Africa <- readRDS('~/G2G_TB/scratch/H3Africa/G2G_Obj_H3Africa.rds')
# G2G_Obj_Tanz <- readRDS('~/G2G_TB/scratch/Tanz/G2G_Obj_Tanz.rds')
# G2G_Obj_TopMed <- readRDS('~/G2G_TB/scratch/Tanz/G2G_Obj_TopMed.rds')
# G2G_Obj_WGS_AFGR <- readRDS('~/G2G_TB/scratch/WGS_AFGR/G2G_Obj_WGS_AFGR.rds')
# G2G_Obj_AFGR_Tanz <- readRDS('~/G2G_TB/scratch/AFGR_Tanz/G2G_Obj_AFGR_Tanz.rds')
G2G_Obj_AFGR_Tanz_ChrX <- readRDS('~/G2G_TB/scratch/AFGR_Tanz_ChrX/G2G_Obj_AFGR_Tanz_ChrX.rds')
G2G_Obj_AFGR_Tanz_ChrX_SIFT <- readRDS('~/G2G_TB/scratch/AFGR_Tanz_ChrX/G2G_Obj_AFGR_Tanz_ChrX_SIFT.rds')

# SummaryStatsTbl(G2G_Obj_AFGR_Tanz)
hits_L3 <- c('Rv3909_Rv3909_p.Val372Ile',"Rv2078_Rv2078_p.Gly74Asp","Rv1516c_Rv1516c_p.Pro37Gln")
hits_L4 <- c("aceAa_Rv1915_p.Asp179Gly","cstA_Rv3063_p.Ser559Arg")

# RunG2G(G2G_Obj_AFGR_Tanz_ChrX,SOFTWARE_DIR,'/home/zmxu/G2G_TB/results/',Ref_Panel = 'AFGR_Tanz_ChrX',tool = 'PLINK',n_cores = 22,n_PC = 3,n_pPC = 0,debug = F,lineage = 'L3',model = NA,X_Chr = T)
# RunG2G(G2G_Obj_AFGR_Tanz_ChrX,SOFTWARE_DIR,'/home/zmxu/G2G_TB/results/',Ref_Panel = 'AFGR_Tanz_ChrX',tool = 'PLINK',n_cores = 22,n_PC = 3,n_pPC = 0,debug = F,lineage = 'L4',model = NA,X_Chr = T)

# RunG2G(G2G_Obj_AFGR_Tanz,SOFTWARE_DIR,'/home/zmxu/G2G_TB/results/',Ref_Panel = 'AFGR_Tanz',tool = 'PLINK',n_cores = 22,n_PC = 3,n_pPC = 0,debug = F,lineage = 'L3',covars_to_incl = c('HIV_Status'),chunks = which(colnames(G2G_Obj_AFGR_Tanz$aa_matrix_filt$L3) %in% hits_L3[1]),model = NA)
# RunG2G(G2G_Obj_AFGR_Tanz,SOFTWARE_DIR,'/home/zmxu/G2G_TB/results/',Ref_Panel = 'AFGR_Tanz',tool = 'PLINK-FIRTH',n_cores = 22,n_PC = 3,n_pPC = 0,debug = F,lineage = 'L4',chunks = which(colnames(G2G_Obj_AFGR_Tanz$aa_matrix_filt$L4) %in% hits_L4))

# RunG2G(G2G_Obj_AFGR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/results/',Ref_Panel = 'AFGR',tool = 'GMMAT-SCORE',n_cores = 22,n_PC = 3,n_pPC = 0,debug = F,lineage = 'L3',chunks = which(colnames(G2G_Obj_AFGR$aa_matrix_filt$L3) %in% hits_L3))
# RunG2G(G2G_Obj_AFGR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/results/',Ref_Panel = 'AFGR',tool = 'GMMAT-SCORE',n_cores = 22,n_PC = 3,n_pPC = 0,debug = F,lineage = 'L4',chunks = which(colnames(G2G_Obj_AFGR$aa_matrix_filt$L4) %in% hits_L4))

# RunG2G(G2G_Obj_AFGR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/results/',Ref_Panel = 'AFGR',tool = 'SAIGE',n_cores = 22,n_PC = 3,n_pPC = 0,debug = F,lineage = 'L3',chunks = which(colnames(G2G_Obj_AFGR$aa_matrix_filt$L3) %in% hits_L3))
# RunG2G(G2G_Obj_AFGR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/results/',Ref_Panel = 'AFGR',tool = 'SAIGE',n_cores = 22,n_PC = 3,n_pPC = 0,debug = F,lineage = 'L4',chunks = which(colnames(G2G_Obj_AFGR$aa_matrix_filt$L4) %in% hits_L4))

# RunG2G(G2G_Obj_WGS_AFGR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/results/',Ref_Panel = 'WGS_AFGR',tool = 'PLINK',n_cores = 22,n_PC = 5,n_pPC = 0,debug = F,lineage = 'L3')
# RunG2G(G2G_Obj_WGS_AFGR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/results/',Ref_Panel = 'WGS_AFGR',tool = 'PLINK',n_cores = 22,n_PC = 5,n_pPC = 0,debug = F,lineage = 'L4')

RunG2GBurden(G2G_Obj_AFGR_Tanz_ChrX,SOFTWARE_DIR,'/home/zmxu/G2G_TB/burden_results/',Ref_Panel = 'AFGR_Tanz_ChrX',tool = 'PLINK',n_PC = 3,n_pPC = 2,n_cores = 20,X_Chr=T,MAC_Thresh = 10,binarize = T,incl_syn = F,gene = 'Rv3909_Rv3909')
RunG2GBurden(G2G_Obj_AFGR_Tanz_ChrX_SIFT,SOFTWARE_DIR,'/home/zmxu/G2G_TB/burden_results_SIFT/',Ref_Panel = 'AFGR_Tanz_ChrX',tool = 'PLINK',n_PC = 3,n_pPC = 2,n_cores = 20,X_Chr=T,MAC_Thresh = 10,binarize = T,incl_syn = F,gene = 'Rv3909_Rv3909')
  
RunG2GBurden(G2G_Obj_AFGR_Tanz_ChrX,SOFTWARE_DIR,'/home/zmxu/G2G_TB/burden_results/',Ref_Panel = 'AFGR_Tanz_ChrX',tool = 'PLINK',n_PC = 3,n_pPC = 2,n_cores = 20,X_Chr=T,MAC_Thresh = 10,binarize = T,incl_syn = T,gene = 'Rv3909_Rv3909')
RunG2GBurden(G2G_Obj_AFGR_Tanz_ChrX_SIFT,SOFTWARE_DIR,'/home/zmxu/G2G_TB/burden_results_SIFT/',Ref_Panel = 'AFGR_Tanz_ChrX',tool = 'PLINK',n_PC = 3,n_pPC = 2,n_cores = 20,X_Chr=T,MAC_Thresh = 10,binarize = T,incl_syn = T,gene = 'Rv3909_Rv3909')

# int_results_ALL <- readRDS('~/G2G_TB/interaction_results/AFGR/PLINK/PC_3_pPC_0_cov_PatientSex-Age/LINEAGE_ALL/results.rds')
# int_results_ALL <- int_results_ALL[sapply(int_results_ALL,nrow) > 0]
# variants_ALL <- sapply(names(int_results_ALL),function(x) strsplit(x=x,split = '.tb_score')[[1]][1])
# ind_ALL <- match(variants_ALL,colnames(G2G_Obj_AFGR$aa_matrix_raw))

# RunInteraction(G2G_Obj_AFGR_Tanz_ChrX,SOFTWARE_DIR,'/home/zmxu/G2G_TB/interaction_results/',Ref_Panel = 'AFGR_Tanz_ChrX',tool = 'PLINK',n_cores = 22,covars_to_incl = c('Patient_Sex','Age','HIV_Status'),n_PC = 5,n_pPC = 0,debug = F,lineage = c('ALL'),by_lineage = T)
# RunInteraction(G2G_Obj_AFGR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/interaction_results/',Ref_Panel = 'AFGR',tool = 'ORDINAL',n_cores = 22,covars_to_incl = c('Patient_Sex','Age','HIV_Status'),n_PC = 3,n_pPC = 0,debug = F,lineage = c('ALL'),by_lineage = T)
# RunInteraction(G2G_Obj_AFGR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/interaction_results/',Ref_Panel = 'AFGR',tool = 'ORDINAL_GxE',n_cores = 22,covars_to_incl = c('Patient_Sex','Age','HIV_Status'),n_PC = 3,n_pPC = 0,debug = F,lineage = c('ALL'),by_lineage = T)

# RunInteraction(G2G_Obj_AFGR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/interaction_results/',Ref_Panel = 'AFGR',tool = 'ORDINAL',n_cores = 22,covars_to_incl = c('Patient_Sex','Age'),n_PC = 3,n_pPC = 3,debug = F,lineage = c('ALL'),chunks = ind_ALL)
# RunInteraction(G2G_Obj_AFGR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/interaction_results/',Ref_Panel = 'AFGR',tool = 'ORDINAL_GxE',n_cores = 22,covars_to_incl = c('Patient_Sex','Age'),n_PC = 3,n_pPC = 0,debug = F,lineage = c('ALL'),chunks = ind_ALL)
# RunInteraction(G2G_Obj_AFGR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/interaction_results/',Ref_Panel = 'AFGR',tool = 'ORDINAL_GxE',n_cores = 22,covars_to_incl = c('Patient_Sex','Age'),n_PC = 3,n_pPC = 3,debug = F,lineage = c('ALL'),chunks = ind_ALL)

# RunG2GPerm(G2G_Obj_AFGR,SOFTWARE_DIR,'/home/zmxu/G2G_TB/perm_results/',Ref_Panel = 'AFGR',tool = 'PLINK',n_cores = 22,n_PC = 3,n_pPC = 0,debug = F,lineage = 'L3',chunks = which(colnames(G2G_Obj_AFGR$aa_matrix_filt$L3) %in% hits_L3[1]))

# L3_Results_AFGR <- ViewResults('~/G2G_TB/results/AFGR/PLINK/PC_3_pPC_0/LINEAGE_L3/results.rds',p_thresh = 5e-8,bonf = F)
# View(do.call(rbind,lapply(1:length(L3_Results_AFGR),function(i) {
#   cur_result <- L3_Results_AFGR[[i]]
#   cbind(cur_result,data.frame(AA = rep(names(L3_Results_AFGR)[i],nrow(cur_result))))
#   })) %>% dplyr::arrange('#CHROM','POS'))
# L4_Results_AFGR <- ViewResults('~/G2G_TB/results/AFGR/PLINK/PC_3_pPC_0/LINEAGE_L4/results.rds',p_thresh = 5e-8,bonf = F)
# View(do.call(rbind,lapply(1:length(L4_Results_AFGR),function(i) {
#   cur_result <- L4_Results_AFGR[[i]]
#   cbind(cur_result,data.frame(AA = rep(names(L4_Results_AFGR)[i],nrow(cur_result))))
# })) %>% dplyr::arrange('#CHROM','POS'))


# L3_Results_WGS_AFGR <- ViewResults('~/G2G_TB/results/WGS_AFGR/PLINK/PC_3_pPC_0/LINEAGE_L3/results.rds',p_thresh = 5e-8,bonf = F)
# View(do.call(rbind,lapply(1:length(L3_Results_WGS_AFGR),function(i) {
#   cur_result <- L3_Results_WGS_AFGR[[i]]
#   cbind(cur_result,data.frame(AA = rep(names(L3_Results_WGS_AFGR)[i],nrow(cur_result))))
#   })) %>% dplyr::arrange('#CHROM','POS'))
# L4_Results_AFGR <- ViewResults('~/G2G_TB/results/AFGR/PLINK/PC_3_pPC_0/LINEAGE_L4/results.rds',p_thresh = 5e-8,bonf = F)
# View(do.call(rbind,lapply(1:length(L4_Results_AFGR),function(i) {
#   cur_result <- L4_Results_AFGR[[i]]
#   cbind(cur_result,data.frame(AA = rep(names(L4_Results_AFGR)[i],nrow(cur_result))))
# })) %>% dplyr::arrange('#CHROM','POS'))


# L3_Results_AFGR_Tanz <- ViewResults('~/G2G_TB/results/AFGR_Tanz/PLINK/PC_3_pPC_0/LINEAGE_L3/results.rds',p_thresh = 5e-8,bonf = F)
# L3_Results_AFGR_Tanz <- do.call(rbind,lapply(1:length(L3_Results_AFGR_Tanz),function(i) {
#   cur_result <- L3_Results_AFGR_Tanz[[i]]
#   cbind(cur_result,data.frame(AA = rep(names(L3_Results_AFGR_Tanz)[i],nrow(cur_result))))
# })) %>% dplyr::arrange('AA')
# 
# L4_Results_AFGR_Tanz  <- ViewResults('~/G2G_TB/results/AFGR_Tanz/PLINK/PC_3_pPC_0/LINEAGE_L4/results.rds',p_thresh = 5e-8,bonf = F)
# L4_Results_AFGR_Tanz <- do.call(rbind,lapply(1:length(L4_Results_AFGR_Tanz),function(i) {
#   cur_result <- L4_Results_AFGR_Tanz[[i]]
#   cbind(cur_result,data.frame(AA = rep(names(L4_Results_AFGR_Tanz)[i],nrow(cur_result))))
# })) %>% dplyr::arrange('AA')

# MAF_L1 <- apply(G2G_Obj_AFGR$aa_matrix_filt$L1,2,function(x) ifelse(length(unique(x)) != 1,min(table(x)),nrow(aa_matrix_raw) - min(table(x))))
# MAF_L2 <- apply(G2G_Obj_AFGR$aa_matrix_filt$L2,2,function(x) ifelse(length(unique(x)) != 1,min(table(x)),nrow(aa_matrix_raw) - min(table(x))))
# MAF_L3 <- apply(G2G_Obj_AFGR$aa_matrix_filt$L3,2,function(x) ifelse(length(unique(x)) != 1,min(table(x)),nrow(aa_matrix_raw) - min(table(x))))
# MAF_L4 <- apply(G2G_Obj_AFGR$aa_matrix_filt$L4,2,function(x) ifelse(length(unique(x)) != 1,min(table(x)),nrow(aa_matrix_raw) - min(table(x))))

# par(mfrow = c(2,2))
# hist(MAC_L1,main = 'Lineage L1',xlab = 'Minor Allele Count')
# hist(MAC_L2,main = 'Lineage L2',xlab = 'Minor Allele Count')
# hist(MAC_L3,main = 'Lineage L3',xlab = 'Minor Allele Count')
# hist(MAC_L4,main = 'Lineage L4',xlab = 'Minor Allele Count')


# pPC1_lineage4 <- sapply(1:ncol(G2G_Obj_AFGR$aa_matrix_filt$L3),function(i) summary(glm(pheno ~ PC1 + PC2 + PC3,data = data.frame(pheno = G2G_Obj_AFGR$aa_matrix_filt$L3[,i],G2G_Obj_AFGR$vir_pPCs$L3)))$coefficients['PC1','Pr(>|t|)'])
# pPC2_lineage4 <- sapply(1:ncol(G2G_Obj_AFGR$aa_matrix_filt$L3),function(i) summary(glm(pheno ~ PC1 + PC2 + PC3,data = data.frame(pheno = G2G_Obj_AFGR$aa_matrix_filt$L3[,i],G2G_Obj_AFGR$vir_pPCs$L3)))$coefficients['PC2','Pr(>|t|)'])
# pPC3_lineage4 <- sapply(1:ncol(G2G_Obj_AFGR$aa_matrix_filt$L3),function(i) summary(glm(pheno ~ PC1 + PC2 + PC3,data = data.frame(pheno = G2G_Obj_AFGR$aa_matrix_filt$L3[,i],G2G_Obj_AFGR$vir_pPCs$L3)))$coefficients['PC3','Pr(>|t|)'])

# par(mfrow = c(2,2))
# plot(1:length(G2G_Obj_AFGR$raw_pPCA$L1$eig),G2G_Obj_AFGR$raw_pPCA$L1$eig / sum(G2G_Obj_AFGR$raw_pPCA$L1$eig),xlab='Number of pPCs',ylab='Variance Explained',main = 'Lineage L1',xlim = c(0,10))
# plot(1:length(G2G_Obj_AFGR$raw_pPCA$L2$eig),G2G_Obj_AFGR$raw_pPCA$L2$eig / sum(G2G_Obj_AFGR$raw_pPCA$L2$eig),xlab='Number of pPCs',ylab='Variance Explained',main = 'Lineage L2',xlim = c(0,10))
# plot(1:length(G2G_Obj_AFGR$raw_pPCA$L3$eig),G2G_Obj_AFGR$raw_pPCA$L3$eig / sum(G2G_Obj_AFGR$raw_pPCA$L3$eig),xlab='Number of pPCs',ylab='Variance Explained',main = 'Lineage L3',xlim = c(0,10))
# plot(1:length(G2G_Obj_AFGR$raw_pPCA$L4$eig),G2G_Obj_AFGR$raw_pPCA$L4$eig / sum(G2G_Obj_AFGR$raw_pPCA$L4$eig),xlab='Number of pPCs',ylab='Variance Explained',main = 'Lineage L4',xlim = c(0,10))


# p_values_min <- sapply(1:60,function(x) min(data.table::fread(cmd = glue::glue("zcat ~/G2G_TB/perm_results/AFGR/PLINK/PC_3_pPC_0/LINEAGE_L3/perm{x}/Rv1516c_Rv1516c_p.Pro37Gln.pval.gz"))$P))

# pPC0_Interaction <- readRDS('~/G2G_TB/interaction_results/AFGR/PLINK/PC_3_pPC_0_cov_PatientSex-Age/LINEAGE_ALL/results.rds')
# pPC0_Interaction <- pPC0_Interaction[sapply(pPC0_Interaction,nrow) > 0]
# pPC3_Interaction <- readRDS('~/G2G_TB/interaction_results/AFGR/PLINK/PC_3_pPC_3_cov_PatientSex-Age/LINEAGE_ALL/results.rds')
# pPC3_Interaction <- pPC3_Interaction[sapply(pPC3_Interaction,nrow) > 0]
# pPC5_Interaction <- readRDS('~/G2G_TB/interaction_results/AFGR/PLINK/PC_3_pPC_5_cov_PatientSex-Age/LINEAGE_ALL/results.rds')
# pPC5_Interaction <- pPC5_Interaction[sapply(pPC5_Interaction,nrow) > 0]

# p_vals <- readRDS('~/G2G_TB/perm_results/AFGR/PLINK/PC_5_pPC_0/LINEAGE_L3/results.rds')
# p_by_variant <- do.call(rbind,p_vals)
# MAC <- apply(G2G_Obj_AFGR$aa_matrix_filt$L3,2,function(x) ifelse(length(unique(x)) != 1,min(table(x)),nrow(AA_Table_filt) - min(table(x))))
# pc_loading <- abs(prcomp(G2G_Obj_AFGR$aa_matrix_filt$L3)$rotation[,1])
# plot(MAC,-log10(apply(p_by_variant,2,median)))
# 
# p_val_Rv3909 <- data.table::fread(cmd = 'zcat ~/G2G_TB/perm_results/AFGR/PLINK/PC_3_pPC_0/LINEAGE_L3/perm3/Rv0004_Rv0004_p.Arg47Trp.pval.gz')
# maf <- data.table::fread('~/G2G_TB/scratch/AFGR/LINEAGE_L3/test.frq')
# 
# d <- ggplot(data.frame(MAF = maf$MAF,P = -log10(p_val_Rv3909$P)), aes(x=MAF, y = P))
# d + geom_hex()

# coef_pPC1 <- sapply(1:ncol(G2G_Obj_AFGR$aa_matrix_filt$L3),function(i) cor(G2G_Obj_AFGR$aa_matrix_filt$L3[,i],G2G_Obj_AFGR$vir_pPCs$L3$PC1,method = 'kendall'))

# p_vals_PC5 <- readRDS('~/G2G_TB/perm_results/AFGR/PLINK/PC_5_pPC_0/LINEAGE_L3/results.rds')
# p_by_variant_PC5 <- do.call(rbind,p_vals_PC5)
# p_values_min <- sapply(1:60,function(x) min(data.table::fread(cmd = glue::glue("zcat ~/G2G_TB/perm_results/AFGR/GMMAT-SCORE/PC_3_pPC_0/LINEAGE_L3/perm{x}/Rv3909_Rv3909_p.Val372Ile.pval.gz"))$PVAL))

# p_vals_PC5_pPC1 <- readRDS('~/G2G_TB/perm_results/AFGR/PLINK/PC_5_pPC_1/LINEAGE_L3/results.rds')
# p_vals_PC5_pPC1 <- lapply(p_vals_PC5_pPC1,function(x){
#   num_vect <- as.numeric(x)
#   names(num_vect) <- names(x)
#   return(num_vect)
# })
# p_by_variant_PC5_pPC1 <- do.call(rbind,p_vals_PC5_pPC1)
# p_values_min <- sapply(1:60,function(x) min(data.table::fread(cmd = glue::glue("zcat ~/G2G_TB/perm_results/AFGR/GMMAT-SCORE/PC_3_pPC_0/LINEAGE_L3/perm{x}/Rv3909_Rv3909_p.Val372Ile.pval.gz"))$PVAL))

# afgr_tanz_dosage <- as(snpStats::read.plink(bed = '/home/zmxu/G2G_TB/scratch/AFGR_Tanz/LINEAGE_L3/TB_DAR_Imputed_G2G.bed',select.snps = 'rs7601374')$genotypes,Class = 'numeric')
# afgr_dosage <- as(snpStats::read.plink(bed = '/home/zmxu/G2G_TB/scratch/AFGR/LINEAGE_L3/TB_DAR_Imputed_G2G.bed',select.snps = 'rs7601374')$genotypes,Class = 'numeric')
# tanz_dosage <- as(snpStats::read.plink(bed = '/home/zmxu/G2G_TB/scratch/Tanz/LINEAGE_L3/TB_DAR_Imputed_G2G.bed',select.snps = '2:160754050:T:C')$genotypes,Class = 'numeric')
# tanz_dosage <- tanz_dosage[match(rownames(afgr_dosage),rownames(tanz_dosage)),]

# H3Africa_dosage <- as(snpStats::read.plink(bed = '/home/zmxu/G2G_TB/scratch/H3Africa/LINEAGE_L3/TB_DAR_Imputed_G2G.bed',select.snps = '2:160754050:T:C')$genotypes,Class = 'numeric')
# H3Africa_dosage <- H3Africa_dosage[match(rownames(afgr_dosage),rownames(H3Africa_dosage)),]

# genotyped_dosage <- as(snpStats::read.plink(bed = '/home/zmxu/G2G_TB/data/Genotyping/TB_DAR_Genotyping/PLINK_QCed_Autosomes_Top/VCF/Fellay_0620.GRCh37.nodup.hwe',select.snps = '2:160754050[b37]T,C')$genotypes,Class = 'numeric')
# genotyped_dosage <- genotyped_dosage[match(rownames(afgr_dosage),rownames(genotyped_dosage)),]

# genotyped_wgs_dosage <- as(snpStats::read.plink(bed = '/home/zmxu/G2G_TB/scratch/WGS_AFGR/LINEAGE_L3/TB_DAR_Imputed_G2G',select.snps = '2:160754050[b37]T,C')$genotypes,Class = 'numeric')

# int_result_L3 <- readRDS('~/G2G_TB/interaction_results/AFGR_Tanz/PLINK/PC_3_pPC_0_cov_PatientSex-Age/LINEAGE_L3/results.rds')
# int_result_L3 <- int_result_L3[sapply(int_result_L3,nrow) > 0]
# int_result_L3 <- do.call(rbind,lapply(1:length(int_result_L3),function(i) {
#   cur_result <- int_result_L3[[i]]
#   cbind(cur_result,data.frame(AA = rep(names(int_result_L3)[i],nrow(cur_result))))
# })) %>% dplyr::arrange('AA')
# 
# int_result_L4 <- readRDS('~/G2G_TB/interaction_results/AFGR_Tanz/PLINK/PC_3_pPC_0_cov_PatientSex-Age/LINEAGE_L4/results.rds')
# int_result_L4 <- int_result_L4[sapply(int_result_L4,nrow) > 0]
# int_result_L4 <- do.call(rbind,lapply(1:length(int_result_L4),function(i) {
#   cur_result <- int_result_L4[[i]]
#   cbind(cur_result,data.frame(AA = rep(names(int_result_L4)[i],nrow(cur_result))))
# })) %>% dplyr::arrange('AA')

# results <- GetResults('~/G2G_TB/results/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_0/LINEAGE_L4/',suffix = '.glm.logistic.hybrid.gz',p_thresh = 5e-8,n_cores = 28,is_interaction = F,is_ordinal = F)
# saveRDS(results,glue::glue("~/G2G_TB/results/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_0/LINEAGE_L4/results.rds"))


Burden_Results_SIFT <- ViewResults('~/G2G_TB/burden_results_SIFT/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_2_Syn/LINEAGE_ALL/results.rds',p_thresh = 5e-8 / length(readRDS('~/G2G_TB/burden_results_SIFT/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_2_Syn/LINEAGE_ALL/results.rds')),bonf = F)
View(do.call(rbind,lapply(1:length(Burden_Results_SIFT),function(i) {
  cur_result <- Burden_Results_SIFT[[i]]
  cbind(cur_result,data.frame(AA = rep(names(Burden_Results_SIFT)[i],nrow(cur_result))))
  })) %>% dplyr::arrange('P'))
