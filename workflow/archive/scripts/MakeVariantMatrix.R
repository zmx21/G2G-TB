library(glue)
library(stringr)
library(ape)
library(parallel)
library(pbmcapply)
library(filematrix)

#Filter variants which only appear in one lineage, and show no variability within the lineage.
RemoveStratVariants <- function(AA_Matrix,Lineage_Df){
  samples_to_keep <- intersect(Lineage_Df$G_NUMBER,rownames(AA_Matrix))
  Lineage_Df <- Lineage_Df[match(samples_to_keep,Lineage_Df$G_NUMBER),]
  #Keep samples where Host data exists
  AA_Matrix_filt <- AA_Matrix[samples_to_keep,]
  #Remove Missing Variants/Genes
  AA_Matrix_filt <- AA_Matrix_filt[,!apply(AA_Matrix_filt,2,function(x) all(is.na(x)))]
  
  #Assign variants to lineages
  strat_table <- lapply(1:ncol(AA_Matrix_filt),function(i) table(ifelse(AA_Matrix_filt[,i]==0,0,1),Lineage_Df$LINEAGE))
  strat_assng <- vector(mode = 'logical',length = length(strat_table))
  for(i in 1:length(strat_table)){
    cur_strat <- strat_table[[i]]
    MAC_by_lineage <- apply(cur_strat,2,function(x) min(x,na.rm = T))
    #If a variant is perfectly stratified by lineage and there are no variability within a lineage, exclude from analysis
    if(sum(MAC_by_lineage > 0) == 0){
      strat_assng[i] <- F
    }else
      strat_assng[i] <- T
  }
  
  return(AA_Matrix_filt[,strat_assng])
}

#Merge AA table of each Mtb sequence into a matrix
AATblToMatrix <- function(AA_Tbl_Files,phylo_snps = NULL,sift_table = NULL,sift_thresh = 0.1,n_cores = 10,missing_matrix){
  if(!is.null(sift_table)){
    sift_excl <- dplyr::filter(sift_table,SIFT_score > sift_thresh)
  }
  #Parse all sample IDs
  all_sample_ids <- sapply(AA_Tbl_Files,function(x) strsplit(x=x,split = '/')[[1]][length(strsplit(x=x,split = '/')[[1]])])
  all_sample_ids <- sapply(all_sample_ids,function(x) gsub(pattern = '.txt',replacement = '',x=x))
  
  #Order missing matrix
  missing_matrix <- missing_matrix[,all_sample_ids]
  
  #Get all non-syn variants for each Mtb sample
  all_variants <- mclapply(AA_Tbl_Files,function(x){
    tbl <- data.table::fread(x)
    if(nrow(tbl) == 0){
      return(data.frame(ID = character(),Genotype = numeric()))
    }
    if(!is.null(phylo_snps)){
      tbl <- dplyr::anti_join(tbl,phylo_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT'))
    }
    if(!is.null(sift_table)){
      tbl <- dplyr::anti_join(tbl,sift_excl,by=c('POS'='POS','REF'='REF','ALT'='ALT'))
    }
    return(data.frame(ID = paste0(tbl$GENE,':',tbl$POS,':',tbl$AA_Change),Genotype = ifelse(tbl$GENOTYPE=='1/1',2,1))) #Homozygotes as 2, heterzygotes (mixed calls) as 1
  },mc.cores = n_cores)
  names(all_variants) <- all_sample_ids
  
  #Construct AA dosage matrix (0,1,or s)
  uniq_variant_ids <- unique(unlist(lapply(all_variants,function(x) x$ID)))
  uniq_variant_genes <- sapply(uniq_variant_ids,function(x) strsplit(x=x,split = ':')[[1]][1])
  uniq_variant_nuc_pos <- as.numeric(sapply(uniq_variant_ids,function(x) strsplit(x=x,split = ':')[[1]][2]))
  uniq_variant_pos <- sapply(uniq_variant_ids,function(x) str_extract(strsplit(x=x,split = 'p\\.')[[1]][2],"[[:digit:]]+"))
  uniq_variant_df <- data.frame(ID = uniq_variant_ids,Gene = uniq_variant_genes, Pos = as.numeric(uniq_variant_pos),Nuc_pos = uniq_variant_nuc_pos)
  uniq_variant_df_sorted <- dplyr::group_by(uniq_variant_df,Gene) %>% dplyr::arrange(Pos,.by_group = T)
  sorted_variant_ids <- uniq_variant_df_sorted$ID
  
  #Initialize AA Matrix
  AA_Matrix <- matrix(0,nrow = length(all_variants),ncol = length(sorted_variant_ids))
  rownames(AA_Matrix) <- names(all_variants)
  colnames(AA_Matrix) <- sorted_variant_ids
  
  #Fill in AA Matrix
  for(i in 1:nrow(AA_Matrix)){
    AA_Matrix[i,all_variants[[i]]$ID] <- all_variants[[i]]$Genotype
    #Set missing positions
    missing_pos <- rownames(missing_matrix)[missing_matrix[,i]==T]
    missing_variants <- dplyr::filter(uniq_variant_df_sorted,Nuc_pos %in% missing_pos)
    AA_Matrix[i,missing_variants$ID] <- NA
  }
  return(AA_Matrix)
}

#Merge Synonymous Nuc table of each Mtb sequence into a matrix
NucSynTblToMatrix <- function(AA_Tbl_Files,phylo_snps = NULL,missing_matrix){
  #Parse all sample IDs
  all_sample_ids <- sapply(AA_Tbl_Files,function(x) strsplit(x=x,split = '/')[[1]][length(strsplit(x=x,split = '/')[[1]])])
  all_sample_ids <- sapply(all_sample_ids,function(x) gsub(pattern = '.txt',replacement = '',x=x))
  
  #Get all non-syn variants for each Mtb sample
  all_variants <- mclapply(AA_Tbl_Files,function(x){
    tbl <- data.table::fread(x)
    if(nrow(tbl) == 0){
      return(data.frame(ID = character(),Genotype = numeric()))
    }
    
    if(!is.null(phylo_snps)){
      tbl <- dplyr::anti_join(tbl,phylo_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT'))
    }
    
    return(data.frame(ID = paste0(tbl$GENE,':',tbl$POS,':',tbl$REF,'>',tbl$ALT),Genotype = ifelse(tbl$GENOTYPE=='1/1',2,1)))
  },mc.cores = n_cores)
  names(all_variants) <- all_sample_ids
  
  #Construct AA dosage matrix (0s or 1s)
  uniq_variant_ids <- unique(unlist(lapply(all_variants,function(x) x$ID)))
  uniq_variant_genes <- sapply(uniq_variant_ids,function(x) strsplit(x=x,split = ':')[[1]][1])
  uniq_variant_pos <- sapply(uniq_variant_ids,function(x) strsplit(x=x,split = ':')[[1]][2])
  uniq_variant_df <- data.frame(ID = uniq_variant_ids,Gene = uniq_variant_genes, Pos = as.numeric(uniq_variant_pos))
  uniq_variant_df_sorted <- dplyr::group_by(uniq_variant_df,Gene) %>% dplyr::arrange(Pos,.by_group = T)
  sorted_variant_ids <- uniq_variant_df_sorted$ID
  
  #Initialize AA Matrix
  AA_Matrix <- matrix(0,nrow = length(all_variants),ncol = length(sorted_variant_ids))
  rownames(AA_Matrix) <- names(all_variants)
  colnames(AA_Matrix) <- sorted_variant_ids
  
  
  #Fill in AA Matrix
  for(i in 1:nrow(AA_Matrix)){
    AA_Matrix[i,all_variants[[i]]$ID] <- all_variants[[i]]$Genotype
    #Set missing positions
    missing_pos <- rownames(missing_matrix)[missing_matrix[,i]==T]
    missing_variants <- dplyr::filter(uniq_variant_df_sorted,Pos %in% missing_pos)
    AA_Matrix[i,missing_variants$ID] <- NA
  }
  return(AA_Matrix)
}

GetGeneBurden <- function(AA_Matrix,AA_Matrix_Syn){
  #Remove SNPs that are lineage markers
  AA_Matrix <- RemoveStratVariants(AA_Matrix,GetParsedMapping())
  
  #Get Variant-Gene mapping for all non-syn variants
  non_syn_variants <- colnames(AA_Matrix)
  non_syn_genes <- sapply(non_syn_variants,function(x) strsplit(x=x,split = ':')[[1]][1])
  non_syn_variant_df <- data.frame(SNP = non_syn_variants,Gene = non_syn_genes,stringsAsFactors = F)
  
  #Get non-syn burden per Gene
  non_syn_genes_uniq <- sort(unique(non_syn_variant_df$Gene))
  non_syn_burden <- matrix(data = NA,nrow = nrow(AA_Matrix),ncol = length(non_syn_genes_uniq))
  rownames(non_syn_burden) <- rownames(AA_Matrix)
  colnames(non_syn_burden) <- non_syn_genes_uniq
  non_syn_burden_homo <- non_syn_burden
  
  for(i in 1:length(non_syn_genes_uniq)){
    cur_gene <- non_syn_genes_uniq[i]
    cur_burden <- apply(AA_Matrix[,dplyr::filter(non_syn_variant_df,Gene == cur_gene)$SNP,drop = F],1,function(x) sum(x != 0,na.rm = F))
    non_syn_burden[,i] <- cur_burden
  }
  
  #Get Variant-Gene mapping for all syn variants, where gene contains syn variants
  syn_variants <- colnames(AA_Matrix_Syn)
  syn_genes <- sapply(syn_variants,function(x) strsplit(x=x,split = ':')[[1]][1])
  syn_variant_df <- data.frame(SNP = syn_variants,Gene = syn_genes) %>% dplyr::filter(Gene %in% non_syn_genes_uniq)
  
  #Get syn burden per Gene
  syn_genes_uniq <- sort(unique(syn_variant_df$Gene))
  syn_burden <- matrix(data = NA,nrow = nrow(AA_Matrix_Syn),ncol = length(syn_genes_uniq))
  rownames(syn_burden) <- rownames(AA_Matrix_Syn)
  colnames(syn_burden) <- syn_genes_uniq
  syn_burden_homo <- syn_burden
  
  for(i in 1:length(syn_genes_uniq)){
    cur_gene <- syn_genes_uniq[i]
    cur_burden <- apply(AA_Matrix_Syn[,dplyr::filter(syn_variant_df,Gene == cur_gene)$SNP,drop = F],1,function(x) sum(x!=0,na.rm = T))
    syn_burden[,i] <- cur_burden

  }
  
  return(list(non_syn_burden=non_syn_burden,syn_burden=syn_burden,non_syn_variant_df=non_syn_variant_df,syn_variant_df=syn_variant_df))
}

args <- commandArgs(trailingOnly = TRUE) 
AA_variant_files <- args[grepl(args,pattern = 'AA_')]
Syn_variant_files <- args[grepl(args,pattern = 'Syn_')]
params <- args[!grepl(args,pattern = 'AA_') & !grepl(args,pattern = 'Syn_')]

missing_mat <- readRDS(params[[1]])
phylo_snps <- params[[2]]
del_tbl <- data.table::fread(params[[3]])
sift_path <- params[[4]]
IS_SIFT <- as.logical(params[[5]])
IS_Burden <- as.logical(params[[6]])
IS_Deletion <- as.logical(params[[7]])
out_path <- params[[8]]
n_cores <- as.numeric(params[[9]])

if(is.na(IS_SIFT) | is.na(IS_Burden) | is.na(IS_Deletion)){
  stop('Please specify options')
}
phylo_snps <- data.table::fread(phylo_snps) %>%
  dplyr::filter(PhylogeneticSNP %in% c('PhyloMarkerANCL1','PhyloMarkerANCL234','PhyloMarkerL1','PhyloMarkerL2','PhyloMarkerL3','PhyloMarkerL4')) %>%
  dplyr::select(POS = Position_ref,REF = anc, ALT = derived)

if(!IS_Burden){
  #Convert AA Tables into a Matrix
  AA_Matrix <- AATblToMatrix(AA_Tbl_Files = AA_variant_files,phylo_snps=phylo_snps,missing_matrix=missing_mat,n_cores = n_cores)
  #Add deleted genes as variants
  if(IS_Deletion){
    del_matrix <- as.matrix(del_tbl[,-c('GNUMBER','LINEAGE','NUMBER_OF_DELETED_GENES')])
    del_matrix[del_matrix==1] <- 2 #Set all deletions as Homo calls
    colnames(del_matrix) <- paste0(colnames(del_matrix),':NA:del')
    rownames(del_matrix) <- del_tbl$GNUMBER
    
    matrix_to_append <- matrix(NA,nrow = nrow(AA_Matrix),ncol = ncol(del_matrix))
    rownames(matrix_to_append) <- rownames(AA_Matrix)
    colnames(matrix_to_append) <- colnames(del_matrix)
    
    matrix_to_append[del_tbl$GNUMBER,] <- del_matrix
    AA_Matrix <- cbind(AA_Matrix,matrix_to_append)
  }
}else if(IS_Burden){
  if(IS_SIFT){
    sift_table <- data.table::fread(cmd = glue::glue("zcat {sift_path}")) %>%
    AA_Matrix_NonSyn <- AATblToMatrix(AA_Tbl_Files = AA_variant_files,phylo_snps=phylo_snps,sift_table = sift_table,missing_matrix=missing_mat,n_cores = n_cores)
  }else{
    AA_Matrix_NonSyn <- AATblToMatrix(AA_Tbl_Files = AA_variant_files,phylo_snps=phylo_snps,missing_matrix=missing_mat,n_cores = n_cores)
  }
  AA_Matrix_Syn <- NucSynTblToMatrix(AA_Tbl_Files = Syn_variant_files,phylo_snps=phylo_snps,missing_matrix = missing_mat,n_cores = n_cores)
  if(IS_Deletion){
    del_matrix <- as.matrix(del_tbl[,-c('GNUMBER','LINEAGE','NUMBER_OF_DELETED_GENES')])
    del_matrix[del_matrix==1] <- 2 #Set all deletions as Homo calls
    
    colnames(del_matrix) <- paste0(colnames(del_matrix),':NA:del')
    rownames(del_matrix) <- del_tbl$GNUMBER
    
    matrix_to_append <- matrix(NA,nrow = nrow(AA_Matrix_NonSyn),ncol = ncol(del_matrix))
    rownames(matrix_to_append) <- rownames(AA_Matrix_NonSyn)
    colnames(matrix_to_append) <- colnames(del_matrix)
    
    matrix_to_append[del_tbl$GNUMBER,] <- del_matrix
    AA_Matrix_NonSyn <- cbind(AA_Matrix_NonSyn,matrix_to_append)
  }
  Gene_Burden <- GetGeneBurden(AA_Matrix = AA_Matrix_NonSyn,AA_Matrix_Syn = AA_Matrix_Syn)
  AA_Matrix <- Gene_Burden
  
}
saveRDS(AA_Matrix,file = out_path)
