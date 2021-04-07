library(glue)
library(stringr)
library(ape)
library(parallel)
library(pbmcapply)
library(filematrix)
DATA_DIR <- '~/G2G_TB/data/'
SOFTWARE_DIR <- '~/G2G_TB/software/'
OUT_DIR <- '~/G2G_TB/scratch/'

Mtb_Data <- paste0(DATA_DIR,"Mtb/Mtb_VCF_files/full/")
Mtb_Data_All_Pos <- paste0(DATA_DIR,"Mtb/Mtb_VCF_files/full_allpos/")
Mtb_Out_AA <- paste0(OUT_DIR,"Mtb_AA_Tbl/")
Mtb_Out_Syn <- paste0(OUT_DIR,"Mtb_Syn_Tbl/")
Mtb_Out_Coverage <- paste0(OUT_DIR,"Mtb_Coverage/")

#Get G_NUMBER, PATIENT_ID and LINEAGE mapping
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

#Convert VCF file to AA list for each sample with non-synonymous variants
WriteAATable <- function(DATA_DIR,SOFTWARE_DIR,OUT_DIR,cur_sample){
  system(glue::glue("mkdir -p {OUT_DIR}"))
  system(glue::glue("mkdir -p {DATA_DIR}/rep_filt/"))
  cur_id <- strsplit(cur_sample,split = '\\.')[[1]][1]
  #Loop through all possible reading frames for each nucleotide variant
  cur_annot_table <- list()
  loop_tbl = T
  cnt <- 0
  #Loop through all frames
  while(loop_tbl){
    #Remove SNPs in repetitive regions 
    system(glue::glue("{SOFTWARE_DIR}bedtools subtract -header -A -a {DATA_DIR}{cur_sample} -b {DATA_DIR}../../Locus_to_exclude_Mtb_Parsed.txt > {DATA_DIR}/rep_filt/{cur_sample}"))
    #Get SNP consequences 
    cur_table <- system(glue::glue("java -jar {SOFTWARE_DIR}snpEff/SnpSift.jar extractFields {DATA_DIR}/rep_filt/{cur_sample} POS REF ALT ANN[{cnt}].EFFECT ANN[{cnt}].GENE EFF[{cnt}].AA > {DATA_DIR}{cur_sample}.tbl"),intern = T)
    cur_table_parsed <- data.table::fread(glue::glue("{DATA_DIR}{cur_sample}.tbl"),sep = '\t',fill = F)
    colnames(cur_table_parsed) <- c('POS','REF','ALT','EFFECT','GENE','AA_Change')
    #End loop if there are no longer have any frames to search on. 
    if(!all(is.na(cur_table_parsed$EFFECT))){
      cnt <- cnt + 1
      cur_annot_table <- c(cur_annot_table,list(cur_table_parsed))
    }else{
      loop_tbl <- F
    }
  }
  #Cleanup
  system(glue::glue("rm {DATA_DIR}{cur_sample}.tbl"))
  if(length(cur_annot_table) == 0){
    return(NA)
  }
  cur_annot_table_full <- do.call(rbind,cur_annot_table) %>% dplyr::filter(EFFECT != 'synonymous_variant' & AA_Change != '')
  data.table::fwrite(cur_annot_table_full,file = glue::glue("{OUT_DIR}{cur_id}.txt"))
}

#Convert VCF file to Nucleotide list for each sample with synonymous variants
WriteNucTableSyn <- function(DATA_DIR,SOFTWARE_DIR,OUT_DIR,cur_sample){
  system(glue::glue("mkdir -p {OUT_DIR}"))
  system(glue::glue("mkdir -p {DATA_DIR}/rep_filt/"))
  cur_id <- strsplit(cur_sample,split = '\\.')[[1]][1]
  #Loop through all possible reading frames for each nucleotide variant
  cur_annot_table <- list()
  loop_tbl = T
  cnt <- 0
  #Loop through all frames
  while(loop_tbl){
    #Remove SNPs in repetitive regions 
    system(glue::glue("{SOFTWARE_DIR}bedtools subtract -header -A -a {DATA_DIR}{cur_sample} -b {DATA_DIR}../../Locus_to_exclude_Mtb_Parsed.txt > {DATA_DIR}/rep_filt/{cur_sample}"))
    #Get SNP consequences 
    cur_table <- system(glue::glue("java -jar {SOFTWARE_DIR}snpEff/SnpSift.jar extractFields {DATA_DIR}/rep_filt/{cur_sample} POS REF ALT ANN[{cnt}].EFFECT ANN[{cnt}].GENE > {DATA_DIR}{cur_sample}.tbl"),intern = T)
    cur_table_parsed <- data.table::fread(glue::glue("{DATA_DIR}{cur_sample}.tbl"),sep = '\t',fill = F)
    colnames(cur_table_parsed) <- c('POS','REF','ALT','EFFECT','GENE')
    #End loop if there are no longer have any frames to search on. 
    if(!all(is.na(cur_table_parsed$EFFECT))){
      cnt <- cnt + 1
      cur_annot_table <- c(cur_annot_table,list(cur_table_parsed))
    }else{
      loop_tbl <- F
    }
  }
  #Cleanup
  system(glue::glue("rm {DATA_DIR}{cur_sample}.tbl"))
  if(length(cur_annot_table) == 0){
    return(NA)
  }
  cur_annot_table_full <- do.call(rbind,cur_annot_table)
  non_syn_variants <- dplyr::filter(cur_annot_table_full,EFFECT %in% c('missense_variant','start_lost','initiator_codon_variant') | grepl(x=EFFECT,pattern = 'stop'))
  syn_variants <- dplyr::filter(cur_annot_table_full,EFFECT == 'synonymous_variant' & !POS %in% non_syn_variants$POS)
  data.table::fwrite(syn_variants,file = glue::glue("{OUT_DIR}{cur_id}.txt"))
}

#Merge AA table of each Mtb sequence into a matrix
AATblToMatrix <- function(AA_Tbl_Files,phylo_snps = NULL,sift_table = NULL,sift_thresh = 0.1,n_cores = 10,cov_thresh = 7,cov_matrix){
  if(!is.null(sift_table)){
    sift_excl <- dplyr::filter(sift_table,SIFT_score > sift_thresh)
  }
  #Parse all sample IDs
  all_sample_ids <- sapply(AA_Tbl_Files,function(x) strsplit(x=x,split = 'Mtb_AA_Tbl/')[[1]][2])
  all_sample_ids <- sapply(all_sample_ids,function(x) gsub(pattern = '.txt',replacement = '',x=x))
  
  #Get all non-syn variants for each Mtb sample
  all_variants <- pbmclapply(AA_Tbl_Files,function(x){
    tbl <- data.table::fread(x)
    if(!is.null(phylo_snps)){
      tbl <- dplyr::anti_join(tbl,phylo_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT'))
    }
    if(!is.null(sift_table)){
      tbl <- dplyr::anti_join(tbl,sift_excl,by=c('POS'='POS','REF'='REF','ALT'='ALT'))
    }
    return(paste0(tbl$GENE,':',tbl$POS,':',tbl$AA_Change))
  },mc.cores = n_cores)
  names(all_variants) <- all_sample_ids
  
  #Construct AA dosage matrix (0s or 1s)
  uniq_variant_ids <- unique(unlist(all_variants))
  uniq_variant_genes <- sapply(uniq_variant_ids,function(x) strsplit(x=x,split = ':')[[1]][1])
  uniq_variant_pos <- sapply(uniq_variant_ids,function(x) str_extract(strsplit(x=x,split = 'p\\.')[[1]][2],"[[:digit:]]+"))
  uniq_variant_df <- data.frame(ID = uniq_variant_ids,Gene = uniq_variant_genes, Pos = as.numeric(uniq_variant_pos))
  uniq_variant_df_sorted <- dplyr::group_by(uniq_variant_df,Gene) %>% dplyr::arrange(Pos,.by_group = T)
  sorted_variant_ids <- uniq_variant_df_sorted$ID
  
  #Initialize AA Matrix
  AA_Matrix <- matrix(0,nrow = length(all_variants),ncol = length(sorted_variant_ids))
  rownames(AA_Matrix) <- names(all_variants)
  colnames(AA_Matrix) <- sorted_variant_ids
  
  #Get coverage for each SNP
  nuc_pos <- as.numeric(sapply(colnames(AA_Matrix),function(x) strsplit(x=x,split = ':')[[1]][2]))
  variant_cov <- t(cov_matrix[nuc_pos,match(rownames(AA_Matrix) ,colnames(cov_matrix))])
  rownames(variant_cov) <- rownames(AA_Matrix) 
  colnames(variant_cov) <- colnames(AA_Matrix)
  
  
  #Fill in AA Matrix
  for(i in 1:nrow(AA_Matrix)){
    AA_Matrix[i,which(cov_matrix[i,] < cov_thresh | is.na(cov_matrix[i,]))] <- NA
    AA_Matrix[i,all_variants[[i]]] <- 1
  }
  return(AA_Matrix)
}

#Merge Synonymous Nuc table of each Mtb sequence into a matrix
NucSynTblToMatrix <- function(AA_Tbl_Files,phylo_snps = NULL,cov_thresh = 7,cov_matrix){
  #Parse all sample IDs
  all_sample_ids <- sapply(AA_Tbl_Files,function(x) strsplit(x=x,split = 'Mtb_Syn_Tbl/')[[1]][2])
  all_sample_ids <- sapply(all_sample_ids,function(x) gsub(pattern = '.txt',replacement = '',x=x))
  
  #Get all non-syn variants for each Mtb sample
  all_variants <- lapply(AA_Tbl_Files,function(x){
    tbl <- data.table::fread(x)
    if(!is.null(phylo_snps)){
      tbl <- dplyr::anti_join(tbl,phylo_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT'))
    }
    
    return(paste0(tbl$GENE,':',tbl$POS,':',tbl$REF,'>',tbl$ALT))
  })
  names(all_variants) <- all_sample_ids
  
  #Construct AA dosage matrix (0s or 1s)
  uniq_variant_ids <- unique(unlist(all_variants))
  uniq_variant_genes <- sapply(uniq_variant_ids,function(x) strsplit(x=x,split = ':')[[1]][1])
  uniq_variant_pos <- sapply(uniq_variant_ids,function(x) strsplit(x=x,split = ':')[[1]][2])
  uniq_variant_df <- data.frame(ID = uniq_variant_ids,Gene = uniq_variant_genes, Pos = as.numeric(uniq_variant_pos))
  uniq_variant_df_sorted <- dplyr::group_by(uniq_variant_df,Gene) %>% dplyr::arrange(Pos,.by_group = T)
  sorted_variant_ids <- uniq_variant_df_sorted$ID

  #Initialize AA Matrix
  AA_Matrix <- matrix(0,nrow = length(all_variants),ncol = length(sorted_variant_ids))
  rownames(AA_Matrix) <- names(all_variants)
  colnames(AA_Matrix) <- sorted_variant_ids
  
  #Get coverage for each SNP
  nuc_pos <- as.numeric(sapply(colnames(AA_Matrix),function(x) strsplit(x=x,split = ':')[[1]][2]))
  variant_cov <- t(cov_matrix[nuc_pos,match(rownames(AA_Matrix) ,colnames(cov_matrix))])
  rownames(variant_cov) <- rownames(AA_Matrix) 
  colnames(variant_cov) <- colnames(AA_Matrix)
  
  
  #Fill in AA Matrix
  for(i in 1:nrow(AA_Matrix)){
    AA_Matrix[i,which(cov_matrix[i,] < cov_thresh | is.na(cov_matrix[i,]))] <- NA
    AA_Matrix[i,all_variants[[i]]] <- 1
  }
  return(AA_Matrix)
}


GetGeneBurden <- function(AA_Matrix,AA_Matrix_Syn,cov_thresh = 7,Gene_Cov_Matrix){
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

  #Get coverage for each gene
  colnames(Gene_Cov_Matrix)[which(grepl(pattern = 'erm',x=colnames(Gene_Cov_Matrix)))] <- as.character(non_syn_genes_uniq[which(grepl(pattern = 'erm',x=non_syn_genes_uniq))])
  Gene_Cov_Matrix_Non_Syn <- Gene_Cov_Matrix[rownames(non_syn_burden),colnames(non_syn_burden)]
  for(i in 1:length(non_syn_genes_uniq)){
    cur_gene <- non_syn_genes_uniq[i]
    cur_burden <- apply(AA_Matrix[,dplyr::filter(non_syn_variant_df,Gene == cur_gene)$SNP,drop = F],1,function(x) sum(x,na.rm = T))
    non_syn_burden[,i] <- cur_burden
    cur_cov <- Gene_Cov_Matrix_Non_Syn[,i]
    non_syn_burden[cur_cov < cov_thresh | is.na(cur_cov),i] <- NA
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
  
  Gene_Cov_Matrix_Syn <- Gene_Cov_Matrix[rownames(non_syn_burden),colnames(non_syn_burden)]
  for(i in 1:length(syn_genes_uniq)){
    cur_gene <- syn_genes_uniq[i]
    cur_burden <- apply(AA_Matrix_Syn[,dplyr::filter(syn_variant_df,Gene == cur_gene)$SNP,drop = F],1,function(x) sum(x,na.rm = T))
    syn_burden[,i] <- cur_burden
    cur_cov <- Gene_Cov_Matrix_Syn[,i]
    syn_burden[cur_cov < cov_thresh | is.na(cur_cov),i] <- NA
  }
  
  return(list(non_syn_burden=non_syn_burden,syn_burden=syn_burden,non_syn_variant_df=non_syn_variant_df,syn_variant_df=syn_variant_df))
}

GetGeneCoverage <- function(SOFTWARE_DIR,Data_Dir,ref_annot,n_cores=60){
  VCF_Files <- dir(Data_Dir)[grepl(pattern = 'all.pos',x=dir(Data_Dir))]
  # system(glue::glue("mkdir -p {Data_Dir}/coverage/"))
  # mclapply(1:length(VCF_Files),function(i){
  #   cur_vcf <- VCF_Files[i]
  #   cur_out <- gsub(x=cur_vcf,pattern = '.vcf.gz','.coverage')
  #   system(glue::glue("{SOFTWARE_DIR}bcftools query -f '%POS [ %SDP]\n' {Data_Dir}/{cur_vcf} > {Data_Dir}/coverage/{cur_out}"))
  # },mc.cores = n_cores)
  
  gff_file <- read.gff(ref_annot)
  gff_file_genes <- gff_file[gff_file$type=='gene',]
  
  gene_df <- data.frame(Gene = sapply(gff_file_genes$attributes,function(x) strsplit(x=strsplit(x=x,split = ';')[[1]][3],split = '=')[[1]][2]),start = gff_file_genes$start, end = gff_file_genes$end,length = gff_file_genes$end - gff_file_genes$start + 1 ,stringsAsFactors = F)

  cov_matrix <- matrix(NA,nrow = length(VCF_Files),ncol = nrow(gene_df))
  sample_names <- sapply(VCF_Files,function(x) strsplit(x=x,split = '\\.')[[1]][1])
  rownames(cov_matrix) <- sample_names
  colnames(cov_matrix) <- gene_df$Gene
  
  for(i in 1:length(VCF_Files)){
    print(i)
    cur_vcf <- VCF_Files[i]
    cur_coverage <- data.table::fread(glue::glue("{Data_Dir}/coverage/{gsub(x=cur_vcf,pattern = '.vcf.gz','.coverage')}"))
    gene_cov <- unlist(mclapply(1:nrow(gene_df),function(k) {
      start <- gene_df$start[k]
      end = gene_df$end[k]
      return(mean(cur_coverage$V2[cur_coverage$V1 >= start & cur_coverage$V1 <= end]))
    },mc.cores = n_cores))
    cov_matrix[i,] <- gene_cov
  }
  return(list(cov_matrix = cov_matrix,gene_df=gene_df))
}

GetSNPCoverage <- function(Data_Dir,n_cores = 20){
  VCF_Files <- dir(Data_Dir)[grepl(pattern = 'all.pos',x=dir(Data_Dir))]
  sample_names <- sapply(VCF_Files,function(x) strsplit(x=x,split = '\\.')[[1]][1])
  
  SNP_Cov <- pbmclapply(VCF_Files,function(x){
    df <- data.table::fread(glue::glue("{Data_Dir}/coverage/{gsub(x=x,pattern = '.vcf.gz','.coverage')}"))
    cov <- rep(NA,length = 4411532)
    cov[df$V1] <- df$V2
    return(cov)
  },mc.cores = n_cores,mc.preschedule = F)
  SNP_Cov_Mat <- do.call(rbind,SNP_Cov)
  rownames(SNP_Cov_Mat) <- sample_names
  return(SNP_Cov_Mat)
}

#Get coverage for each gene
# Gene_Result <- GetGeneCoverage(SOFTWARE_DIR,Mtb_Data_All_Pos,'~/G2G_TB/data/Mtb/Mtb_VCF_files/GCF_000195955.2.gff')
# Gene_Cov_Matrix <- Gene_Resul$cov_matrix
# saveRDS(Gene_Cov_Matrix,glue::glue("{Mtb_Out_Coverage}/Mtb_Gene_Coverage.rds"))
# Gene_Length <- Gene_Result$gene_df
# data.table::fwrite(Gene_Length,glue::glue("{Mtb_Out_Coverage}/Mtb_Gene_Length.txt"),sep = ' ',col.names = T,row.names = F)
# Gene_Cov_Matrix <- readRDS(glue::glue("{Mtb_Out_Coverage}/Mtb_Gene_Coverage.rds"))

# SNP_Cov_Matrix <- GetSNPCoverage(Mtb_Data_All_Pos)
# fm = fm.create.from.matrix(t(SNP_Cov_Matrix),filenamebase = glue::glue("{Mtb_Out_Coverage}/Mtb_SNP_Coverage"))
# SNP_Cov_Matrix <- fm.open(glue::glue("{Mtb_Out_Coverage}/Mtb_SNP_Coverage"))


#Get VCF Files and write into AA tables
# VCF_Files <- dir(Mtb_Data)[grepl(pattern = 'var.homo',x=dir(Mtb_Data))]
# pbmclapply(VCF_Files,function(x) WriteNucTableSyn(Mtb_Data,SOFTWARE_DIR,Mtb_Out_Syn,x),mc.cores = 20)
# pbmclapply(VCF_Files,function(x) WriteAATable(Mtb_Data,SOFTWARE_DIR,Mtb_Out_AA,x),mc.cores = 20)

#Convert AA Tables into a Matrix
# AA_Tbl_Files <- dir(Mtb_Out_AA)[grepl(pattern = '.txt',x=dir(Mtb_Out_AA))]
# AA_Matrix <- AATblToMatrix(paste0(Mtb_Out_AA,AA_Tbl_Files),cov_matrix = SNP_Cov_Matrix)
# saveRDS(AA_Matrix,file = glue::glue("{Mtb_Out_AA}AA_Matrix.rds"))
# AA_Matrix <- readRDS(glue::glue("{Mtb_Out_AA}AA_Matrix.rds"))

#Convert Syn nucleotide tables to matrix
# AA_Tbl_Files_Syn <- dir(Mtb_Out_Syn)[grepl(pattern = '.txt',x=dir(Mtb_Out_Syn))]
# AA_Matrix_Syn <- NucSynTblToMatrix(paste0(Mtb_Out_Syn,AA_Tbl_Files_Syn),cov_matrix = SNP_Cov_Matrix)
# saveRDS(AA_Matrix_Syn,file = glue::glue("{Mtb_Out_Syn}AA_Matrix.rds"))
# AA_Matrix_Syn <- readRDS(glue::glue("{Mtb_Out_Syn}AA_Matrix.rds"))

#Construct Burden Table for each gene
# phylo_snps <- data.table::fread(glue::glue("{DATA_DIR}/pheno/Mtb/PositionsPhylogeneticSNPs_20171004.txt")) %>%
#   dplyr::filter(PhylogeneticSNP %in% c('PhyloMarkerANCL1','PhyloMarkerANCL234','PhyloMarkerL1','PhyloMarkerL2','PhyloMarkerL3','PhyloMarkerL4')) %>%
#   dplyr::select(POS = Position_ref,REF = anc, ALT = derived)
# 
# AA_Matrix_No_Phylo <- AATblToMatrix(paste0(Mtb_Out_AA,AA_Tbl_Files),phylo_snps=phylo_snps,cov_matrix = SNP_Cov_Matrix)
# AA_Matrix_Syn_No_Phylo <- NucSynTblToMatrix(paste0(Mtb_Out_Syn,AA_Tbl_Files_Syn),phylo_snps=phylo_snps,cov_matrix = SNP_Cov_Matrix)
# Gene_Burden <- GetGeneBurden(AA_Matrix = AA_Matrix_No_Phylo,AA_Matrix_Syn = AA_Matrix_Syn_No_Phylo,Gene_Cov_Matrix = Gene_Cov_Matrix)
# saveRDS(Gene_Burden,file = glue::glue("{Mtb_Out_Syn}Gene_Burden.rds"))

# sift_table <- data.table::fread(cmd = glue::glue("zcat {DATA_DIR}/Mtb/SIFT/GCA_000195955.2.22/Chromosome.gz")) %>%
#   dplyr::select(POS = `#Position`,REF = Ref_allele,ALT = New_allele,SIFT_score)
# AA_Matrix_No_Phylo_SIFT <- AATblToMatrix(paste0(Mtb_Out_AA,AA_Tbl_Files),phylo_snps=phylo_snps,sift_table = sift_table,cov_matrix = SNP_Cov_Matrix)
# Gene_Burden_SIFT <- GetGeneBurden(AA_Matrix = AA_Matrix_No_Phylo_SIFT,AA_Matrix_Syn = AA_Matrix_Syn_No_Phylo,Gene_Cov_Matrix = Gene_Cov_Matrix)
# saveRDS(Gene_Burden_SIFT,file = glue::glue("{Mtb_Out_Syn}Gene_Burden_SIFT.rds"))

