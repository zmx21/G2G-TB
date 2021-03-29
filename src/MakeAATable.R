library(glue)
library(stringr)
library(ape)
library(parallel)

DATA_DIR <- '~/G2G_TB/data/'
SOFTWARE_DIR <- '~/G2G_TB/software/'
OUT_DIR <- '~/G2G_TB/scratch/'

Mtb_Data <- paste0(DATA_DIR,"Mtb/Mtb_VCF_files/full/")
Mtb_Data_All_Pos <- paste0(DATA_DIR,"Mtb/Mtb_VCF_files/full_allpos/")
Mtb_Out_AA <- paste0(OUT_DIR,"Mtb_AA_Tbl/")
Mtb_Out_Syn <- paste0(OUT_DIR,"Mtb_Syn_Tbl/")
Mtb_Out_Coverage <- paste0(OUT_DIR,"Mtb_Coverage/")

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
AATblToMatrix <- function(AA_Tbl_Files,phylo_snps = NULL,sift_table = NULL,sift_thresh = 0.1,n_cores = 10){
  sift_excl <- dplyr::filter(sift_table,SIFT_score > sift_thresh)
  #Parse all sample IDs
  all_sample_ids <- sapply(AA_Tbl_Files,function(x) strsplit(x=x,split = 'Mtb_AA_Tbl/')[[1]][2])
  all_sample_ids <- sapply(all_sample_ids,function(x) gsub(pattern = '.txt',replacement = '',x=x))
  
  #Get all non-syn variants for each Mtb sample
  all_variants <- pbmclapply(AA_Tbl_Files,function(x){
    tbl <- data.table::fread(x)
    if(!is.null(phylo_snps)){
      tbl <- dplyr::anti_join(tbl,phylo_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT'))
    }
    if(!is.null(sift_thresh)){
      tbl <- dplyr::anti_join(tbl,sift_excl,by=c('POS'='POS','REF'='REF','ALT'='ALT'))
    }
    return(paste0(tbl$GENE,'_',tbl$AA_Change))
  },mc.cores = n_cores)
  names(all_variants) <- all_sample_ids
  
  #Construct AA dosage matrix (0s or 1s)
  uniq_variant_ids <- unique(unlist(all_variants))
  uniq_variant_genes <- sapply(uniq_variant_ids,function(x) strsplit(x=x,split = 'p\\.')[[1]][1])
  uniq_variant_pos <- sapply(uniq_variant_ids,function(x) str_extract(strsplit(x=x,split = 'p\\.')[[1]][2],"[[:digit:]]+"))
  uniq_variant_df <- data.frame(ID = uniq_variant_ids,Gene = uniq_variant_genes, Pos = as.numeric(uniq_variant_pos))
  uniq_variant_df_sorted <- dplyr::group_by(uniq_variant_df,Gene) %>% dplyr::arrange(Pos,.by_group = T)
  sorted_variant_ids <- uniq_variant_df_sorted$ID
  
  AA_Matrix <- matrix(0,nrow = length(all_variants),ncol = length(sorted_variant_ids))
  rownames(AA_Matrix) <- names(all_variants)
  colnames(AA_Matrix) <- sorted_variant_ids
  #Fill in AA Matrix
  for(i in 1:nrow(AA_Matrix)){
    AA_Matrix[i,all_variants[[i]]] <- 1
  }
  return(AA_Matrix)
}

#Merge Synonymous Nuc table of each Mtb sequence into a matrix
NucSynTblToMatrix <- function(AA_Tbl_Files,phylo_snps){
  #Parse all sample IDs
  all_sample_ids <- sapply(AA_Tbl_Files,function(x) strsplit(x=x,split = 'Mtb_Syn_Tbl/')[[1]][2])
  all_sample_ids <- sapply(all_sample_ids,function(x) gsub(pattern = '.txt',replacement = '',x=x))
  
  #Get all non-syn variants for each Mtb sample
  all_variants <- lapply(AA_Tbl_Files,function(x){
    tbl <- data.table::fread(x)
    if(!is.null(phylo_snps)){
      tbl <- dplyr::anti_join(tbl,phylo_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT'))
    }
    
    return(paste0(tbl$GENE,'_',tbl$POS,':',tbl$REF,'>',tbl$ALT))
  })
  names(all_variants) <- all_sample_ids
  
  #Construct AA dosage matrix (0s or 1s)
  uniq_variant_ids <- unique(unlist(all_variants))
  uniq_variant_genes <- sapply(uniq_variant_ids,function(x) paste0(strsplit(x=x,split = '_')[[1]][1:2],collapse = '_'))
  uniq_variant_pos <- sapply(uniq_variant_ids,function(x) str_extract(strsplit(x=x,split = '_')[[1]][3],"[[:digit:]]+"))
  uniq_variant_df <- data.frame(ID = uniq_variant_ids,Gene = uniq_variant_genes, Pos = as.numeric(uniq_variant_pos))
  uniq_variant_df_sorted <- dplyr::group_by(uniq_variant_df,Gene) %>% dplyr::arrange(Pos,.by_group = T)
  sorted_variant_ids <- uniq_variant_df_sorted$ID
  
  AA_Matrix <- matrix(0,nrow = length(all_variants),ncol = length(sorted_variant_ids))
  rownames(AA_Matrix) <- names(all_variants)
  colnames(AA_Matrix) <- sorted_variant_ids
  #Fill in AA Matrix
  for(i in 1:nrow(AA_Matrix)){
    AA_Matrix[i,all_variants[[i]]] <- 1
  }
  return(AA_Matrix)
}


GetGeneBurden <- function(AA_Matrix,AA_Matrix_Syn){
  #Get Variant-Gene mapping for all non-syn variants
  non_syn_variants <- colnames(AA_Matrix)
  non_syn_genes <- sapply(non_syn_variants,function(x) paste0(strsplit(x=x,split = '_')[[1]][1:2],collapse = '_'))
  non_syn_variant_df <- data.frame(SNP = non_syn_variants,Gene = non_syn_genes)
  
  #Get non-syn burden per Gene
  non_syn_genes_uniq <- sort(unique(non_syn_variant_df$Gene))
  non_syn_burden <- matrix(data = NA,nrow = nrow(AA_Matrix),ncol = length(non_syn_genes_uniq))
  rownames(non_syn_burden) <- rownames(AA_Matrix)
  colnames(non_syn_burden) <- non_syn_genes_uniq
  for(i in 1:length(non_syn_genes_uniq)){
    cur_gene <- non_syn_genes_uniq[i]
    cur_burden <- apply(AA_Matrix[,dplyr::filter(non_syn_variant_df,Gene == cur_gene)$SNP,drop = F],1,sum)
    non_syn_burden[,i] <- cur_burden
  }
  
  #Get Variant-Gene mapping for all syn variants, where gene contains syn variants
  syn_variants <- colnames(AA_Matrix_Syn)
  syn_genes <- sapply(syn_variants,function(x) paste0(strsplit(x=x,split = '_')[[1]][1:2],collapse = '_'))
  syn_variant_df <- data.frame(SNP = syn_variants,Gene = syn_genes) %>% dplyr::filter(Gene %in% non_syn_genes_uniq)
  
  #Get syn burden per Gene
  syn_genes_uniq <- sort(unique(syn_variant_df$Gene))
  syn_burden <- matrix(data = NA,nrow = nrow(AA_Matrix_Syn),ncol = length(syn_genes_uniq))
  rownames(syn_burden) <- rownames(AA_Matrix_Syn)
  colnames(syn_burden) <- syn_genes_uniq
  
  for(i in 1:length(syn_genes_uniq)){
    cur_gene <- syn_genes_uniq[i]
    cur_burden <- apply(AA_Matrix_Syn[,dplyr::filter(syn_variant_df,Gene == cur_gene)$SNP,drop = F],1,sum)
    syn_burden[,i] <- cur_burden
  }
  
  return(list(non_syn_burden=non_syn_burden,syn_burden=syn_burden,non_syn_variant_df=non_syn_variant_df,syn_variant_df=syn_variant_df))
}

GetCoverage <- function(SOFTWARE_DIR,Data_Dir,ref_annot,n_cores=60){
  VCF_Files <- dir(Data_Dir)[grepl(pattern = 'all.pos',x=dir(Data_Dir))]
  # system(glue::glue("mkdir -p {Data_Dir}/coverage/"))
  # mclapply(1:length(VCF_Files),function(i){
  #   cur_vcf <- VCF_Files[i]
  #   cur_out <- gsub(x=cur_vcf,pattern = '.vcf.gz','.coverage')
  #   system(glue::glue("{SOFTWARE_DIR}bcftools query -f '%POS [ %SDP]\n' {Data_Dir}/{cur_vcf} > {Data_Dir}/coverage/{cur_out}"))
  # },mc.cores = n_cores)
  
  gff_file <- read.gff(ref_annot)
  gff_file_genes <- gff_file[gff_file$type=='gene',]
  
  gene_df <- data.frame(Gene = sapply(gff_file_genes$attributes,function(x) strsplit(x=strsplit(x=x,split = ';')[[1]][3],split = '=')[[1]][2]),start = gff_file_genes$start, end = gff_file_genes$end,stringsAsFactors = F)
  
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
  return(cov_matrix)
}
#Get coverage for each gene
Cov_Matrix <- GetCoverage(SOFTWARE_DIR,Mtb_Data_All_Pos,'~/G2G_TB/data/Mtb/Mtb_VCF_files/GCF_000195955.2.gff')
saveRDS(Cov_Matrix,glue::glue("{Mtb_Out_Coverage}/Mtb_Coverage.rds"))

#Get VCF Files and write into AA tables
# VCF_Files <- dir(Mtb_Data)[grepl(pattern = 'var.homo',x=dir(Mtb_Data))]
# pbmclapply(VCF_Files,function(x) WriteNucTableSyn(Mtb_Data,SOFTWARE_DIR,Mtb_Out_Syn,x),mc.cores = 20)
# pbmclapply(VCF_Files,function(x) WriteAATable(Mtb_Data,SOFTWARE_DIR,Mtb_Out_AA,x),mc.cores = 20)

#Convert AA Tables into a Matrix
# AA_Tbl_Files <- dir(Mtb_Out_AA)[grepl(pattern = '.txt',x=dir(Mtb_Out_AA))]
# AA_Matrix <- AATblToMatrix(paste0(Mtb_Out_AA,AA_Tbl_Files))

# saveRDS(AA_Matrix,file = glue::glue("{Mtb_Out_AA}AA_Matrix.rds"))


#Convert Syn nucleotide tables to matrix
# AA_Tbl_Files_Syn <- dir(Mtb_Out_Syn)[grepl(pattern = '.txt',x=dir(Mtb_Out_Syn))]
# AA_Matrix_Syn <- NucSynTblToMatrix(paste0(Mtb_Out_Syn,AA_Tbl_Files_Syn))
# saveRDS(AA_Matrix_Syn,file = glue::glue("{Mtb_Out_Syn}AA_Matrix.rds"))


#Construct Burden Table for each gene
# phylo_snps <- data.table::fread(glue::glue("{DATA_DIR}/pheno/Mtb/PositionsPhylogeneticSNPs_20171004.txt")) %>% 
#   dplyr::filter(PhylogeneticSNP %in% c('PhyloMarkerANCL1','PhyloMarkerANCL234','PhyloMarkerL1','PhyloMarkerL2','PhyloMarkerL3','PhyloMarkerL4')) %>%
#   dplyr::select(POS = Position_ref,REF = anc, ALT = derived)
# AA_Matrix_No_Phylo <- AATblToMatrix(paste0(Mtb_Out_AA,AA_Tbl_Files),phylo_snps)
# AA_Matrix_Syn_No_Phylo <- NucSynTblToMatrix(paste0(Mtb_Out_Syn,AA_Tbl_Files_Syn),phylo_snps)
# Gene_Burden <- GetGeneBurden(AA_Matrix_No_Phylo,AA_Matrix_Syn_No_Phylo)
# saveRDS(Gene_Burden,file = glue::glue("{Mtb_Out_Syn}Gene_Burden.rds"))

# sift_table <- data.table::fread(cmd = glue::glue("zcat {DATA_DIR}/Mtb/SIFT/GCA_000195955.2.22/Chromosome.gz")) %>% 
#   dplyr::select(POS = `#Position`,REF = Ref_allele,ALT = New_allele,SIFT_score)
# AA_Matrix_No_Phylo_SIFT <- AATblToMatrix(paste0(Mtb_Out_AA,AA_Tbl_Files),phylo_snps,sift_table)
# Gene_Burden_SIFT <- GetGeneBurden(AA_Matrix_No_Phylo_SIFT,AA_Matrix_Syn_No_Phylo)
# saveRDS(Gene_Burden_SIFT,file = glue::glue("{Mtb_Out_Syn}Gene_Burden_SIFT.rds"))

