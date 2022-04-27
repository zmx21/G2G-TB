library(dplyr)
library(ape)
library(phylobase)
library(adephylo)

#Minor allele count of a TB variant
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
#Phylogenetic PC 
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

SetUpG2G <- function(OUT_DIR,Metadata_Path,Host_PC_Path,Var_Tbl_Path,Phylo_Tree_Path,Mtb_Nuc_Out,Host_MAF,Pathogen_MAC_pPCA,Pathogen_MAC_AA_Lineage,Pathogen_MAC_AA,BFILE_Path,Host_Files,n_cores){

  #Get metadata file
  parsed_mapping <- data.table::fread(Metadata_Path) %>% dplyr::select(PATIENT_ID,G_NUMBER,LINEAGE=Lineage)
  
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
  host_IDs <- data.table::fread(glue::glue("{BFILE_Path}.fam"))$V2
  host_IDs_simple <- sapply(host_IDs,function(x) gsub(x=x,pattern = 'WGS_Fellay.',replacement = ''))
  host_IDs_simple <- sapply(host_IDs_simple,function(x) gsub(x=x,pattern = 'Batch1_',replacement = ''))
  host_IDs_simple <- sapply(host_IDs_simple,function(x) gsub(x=x,pattern = 'Batch2_',replacement = ''))
  
  #Initialize variables to store
  both_IDs_to_keep <- vector(mode = 'list',length = length(vir_pPCs)); names(both_IDs_to_keep) <- names(vir_pPCs)
  host_PCs <- vector(mode = 'list',length = length(vir_pPCs)); names(host_PCs) <- names(vir_pPCs)
  covars <- vector(mode = 'list',length = length(vir_pPCs)); names(covars) <- names(vir_pPCs)
  
  #Get Host PCs
  host_PCs_raw <- data.table::fread(Host_PC_Path)
  
  #Get Covars
  covars_discrete <- data.table::fread(glue::glue("{Host_Files}/covars_discrete"),sep = '\t',na.strings = c("",'unknown')) 
  covars_numeric <- data.table::fread(glue::glue("{Host_Files}/covars_numeric"),sep = '\t',na.strings = c("",'unknown')) 
  covars_raw <- dplyr::inner_join(covars_discrete,covars_numeric)
  
  #Loop through lineages, store PLINK files, covars, and host PCs
  for(i in 1:length(vir_pPCs)){
    system(glue::glue('mkdir -p {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/'))
    pathogen_IDs_to_keep <- dplyr::filter(parsed_mapping,G_NUMBER %in% vir_pPCs[[i]]$G_NUMBER)
    both_IDs_to_keep[[i]] <- dplyr::filter(pathogen_IDs_to_keep,PATIENT_ID %in% host_IDs_simple)
    both_IDs_to_keep[[i]]$FAM_ID <- host_IDs[match(both_IDs_to_keep[[i]]$PATIENT_ID,host_IDs_simple)]
    data.table::fwrite(data.frame(FID=both_IDs_to_keep[[i]]$FAM_ID,IID=both_IDs_to_keep[[i]]$FAM_ID),
                       file = glue::glue("{OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt"),sep = '\t',quote = F)
    
    #Write out PLINK files with samples to keep
    system(
      glue::glue(
        "~/Software/plink2 --bfile {BFILE_Path} --keep {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --indiv-sort f {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --make-bed --out {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_G2G_Raw"
      )
    )
    
    #Apply MAF Filter
    system(
      glue::glue(
        "~/Software/plink2 --bfile {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_G2G_Raw --maf {Host_MAF} --indiv-sort f {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --make-bed --out {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_G2G"
      )
    )
    
    #Write out HLA PLINK files with samples to keep
    system(
      glue::glue(
        "~/Software/plink2 --bfile {Host_Files}/TB_DAR_HLA_Alleles --keep {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --indiv-sort f {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --make-bed --out {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_HLA_Alleles_G2G"
      )
    )
    system(
      glue::glue(
        "~/Software/plink2 --bfile {Host_Files}/TB_DAR_HLA_AA --keep {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --indiv-sort f {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/G2G_Samples.txt --make-bed --out {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_HLA_AA_G2G"
      )
    )
    
    system(glue::glue("rm {OUT_DIR}/LINEAGE_{names(vir_pPCs)[i]}/TB_DAR_G2G_Raw*"))
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
  aa_matrix_full <- aa_matrix[both_IDs_to_keep[['ALL']]$G_NUMBER,,drop=FALSE]
  aa_matrix_full<- aa_matrix_full[,GetMAC(aa_matrix_full) > Pathogen_MAC_AA,drop=FALSE]
  
  
  #Filter AA Matrix, decide for each variant whether to do stratified or stratified analysis
  aa_matrix_filt <- FilterAAMatrix(aa_matrix_full,both_IDs_to_keep,MAC_Thresh = Pathogen_MAC_AA_Lineage) 
  
  #Path to VCF file (for each lineage)
  host_path <- glue::glue("{OUT_DIR}/LINEAGE_{c(names(vir_pPCs),'ALL')}/TB_DAR_G2G")
  names(host_path) <- c(names(vir_pPCs),'ALL')
  
  #Load TB Score
  tb_score_df <- data.table::fread(glue::glue("{Host_Files}/TB_score"))
  tb_score_by_lineage <- lapply(both_IDs_to_keep,function(x) dplyr::left_join(x %>% dplyr::select(IID=FAM_ID),tb_score_df,c('IID'='IID')) %>% dplyr::relocate(FID,IID,TB_score))
  
  #Load X-Ray Score
  xray_score_df <- data.table::fread(glue::glue("{Host_Files}/Xray_score"))
  xray_score_by_lineage <- lapply(both_IDs_to_keep,function(x) dplyr::left_join(x %>% dplyr::select(IID=FAM_ID),xray_score_df,c('IID'='IID')) %>% dplyr::relocate(FID,IID,Xray_score))
  
  return(list(host_PCs=host_PCs,
              vir_pPCs=vir_pPCs,
              covars = covars,
              tb_score=tb_score_by_lineage,
              xray_score = xray_score_by_lineage,
              aa_matrix_filt=aa_matrix_filt,
              aa_matrix_full=aa_matrix_full,
              nuc_matrix = nuc_matrix_trans,
              both_IDs_to_keep=both_IDs_to_keep,
              host_path = host_path,
              raw_pPCA = vir_pPCA))
}


args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- args[[1]]
Metadata_Path <- args[[2]]
Host_PC_Path <- args[[3]]
Var_Tbl_Path <- args[[4]]
Phylo_Tree_Path <- args[[5]]
Mtb_Nuc_Out <- args[[6]]
Host_MAF <- as.numeric(args[[7]])
Pathogen_MAC_pPCA <- as.numeric(args[[8]])
Pathogen_MAC_AA_Lineage <- as.numeric(args[[9]])
Pathogen_MAC_AA <- as.numeric(args[[10]])
BFILE_Path <- gsub(args[[11]],pattern = '.bed',replacement = '')
Host_Files <- args[[12]]
n_cores <- as.numeric(args[[13]])

# OUT_DIR <- '../../results/Burden_False_SIFT_False_Del_True_HomoOnly_True/'
# Metadata_Path <- '../../data/pheno/metadata_Sinergia_final_dataset_human_bac_genome_available_QCed.txt'
# Host_PC_Path <- '../../scratch/Host/TB_DAR_GWAS_PCA.eigenvec'
# Var_Tbl_Path <- '../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True/Mtb_Var_Tbl.rds'
# Phylo_Tree_Path <- '../../data/Mtb/RAxML_bestTree.Sinergia_final_dataset_human_bac_genome_available_rerooted.nwk'
# Mtb_Nuc_Out <- '../../data/Mtb/merged/merged.var.homo.SNPs.vcf.dosage'
# Host_MAF <- 0.05
# Pathogen_MAC_pPCA <- 10
# Pathogen_MAC_AA_Lineage <- 10
# Pathogen_MAC_AA <- 10
# BFILE_Path <- '../../data/Genotyping_WGS/TBDAR.WGS.Imputed.GWASReady'
# Host_Files <- '../../scratch/Host/'
# n_cores <- 1

G2G_Obj <- SetUpG2G(OUT_DIR=OUT_DIR,Metadata_Path=Metadata_Path,Host_PC_Path=Host_PC_Path,Var_Tbl_Path=Var_Tbl_Path,Phylo_Tree_Path=Phylo_Tree_Path,Mtb_Nuc_Out=Mtb_Nuc_Out,Host_MAF=Host_MAF,Pathogen_MAC_pPCA=Pathogen_MAC_pPCA,Pathogen_MAC_AA_Lineage=Pathogen_MAC_AA_Lineage,Pathogen_MAC_AA=Pathogen_MAC_AA,BFILE_Path=BFILE_Path,Host_Files=Host_Files,n_cores=n_cores)
saveRDS(G2G_Obj,glue::glue("{OUT_DIR}/G2G_Obj.rds"))
