library(ape)
library(parallel)
GetGeneCoverage <- function(VCF_Files,ref_annot,n_cores){
  pbmclapply(1:length(VCF_Files),function(i){
    cur_vcf <- VCF_Files[i]
    cur_out <- gsub(x=cur_vcf,pattern = '.vcf.gz','.coverage')
    system(glue::glue("bcftools query -f '%POS [ %SDP]\n' {cur_vcf} > {cur_out}"))
  },mc.cores = n_cores)
  
  gff_file <- read.gff(ref_annot)
  gff_file_genes <- gff_file[gff_file$type=='gene',]
  
  gene_df <- data.frame(Gene = sapply(gff_file_genes$attributes,function(x) strsplit(x=strsplit(x=x,split = ';')[[1]][3],split = '=')[[1]][2]),start = gff_file_genes$start, end = gff_file_genes$end,length = gff_file_genes$end - gff_file_genes$start + 1 ,stringsAsFactors = F)
  
  cov_matrix <- matrix(NA,nrow = length(VCF_Files),ncol = nrow(gene_df))
  sample_names <- sapply(VCF_Files,function(x) strsplit(x=strsplit(x=x,split = '/')[[1]][length(strsplit(x=x,split = '/')[[1]])],split = '\\.')[[1]][1])
  rownames(cov_matrix) <- sample_names
  colnames(cov_matrix) <- gene_df$Gene
  
  for(i in 1:length(VCF_Files)){
    print(i)
    cur_vcf <- VCF_Files[i]
    cur_out <- gsub(x=cur_vcf,pattern = '.vcf.gz','.coverage')
    
    cur_coverage <- data.table::fread(cur_out)
    gene_cov <- unlist(pbmclapply(1:nrow(gene_df),function(k) {
      start <- gene_df$start[k]
      end = gene_df$end[k]
      return(mean(cur_coverage$V2[cur_coverage$V1 >= start & cur_coverage$V1 <= end]))
    },mc.cores = n_cores))
    cov_matrix[i,] <- gene_cov
    system(glue::glue("rm {cur_out}"))
  }
  return(list(cov_matrix = cov_matrix,gene_df=gene_df))
}

#Get positions for each samples which contain missing call
GetMissingMatrix <- function(missing_files){
  #Get data frame for each sample, containing missing positions
  missing_tbl <- lapply(missing_files,function(x) data.table::fread(x))
  names(missing_tbl) <- sapply(missing_files,function(x) {
    split_path <- strsplit(x=x,split = '/')[[1]]
    file_name <- split_path[length(split_path)]
    ID <- strsplit(file_name,split = '\\.')[[1]][1]
    return(ID)
  })
  
  uniq_pos <- sort(unique(do.call(rbind, missing_tbl)$V1))
  missing_mat <- matrix(F,nrow = length(uniq_pos),ncol = length(missing_tbl))
  rownames(missing_mat) <- as.character(uniq_pos)
  colnames(missing_mat) <- names(missing_tbl)
  
  for(i in 1:ncol(missing_mat)){
    missing_mat[as.character(missing_tbl[[i]]$V1),i] <- T
  }
  return(missing_mat)
}
args <- commandArgs(trailingOnly = TRUE) 
vcf_files <- args[grepl(args,pattern = '.vcf.gz')]
missing_files <- args[grepl(args,pattern = '.missing')]
params <- args[!grepl(args,pattern = '.vcf.gz') & !grepl(args,pattern = '.missing')]
ref_annot <- params[[1]]
out_dir <- params[[2]]
n_cores <- as.numeric(params[[3]])
system(glue::glue("mkdir -p {out_dir}"))

#Get average coverage for each gene
# Gene_Result <- GetGeneCoverage(vcf_files,ref_annot,n_cores)
# Gene_Cov_Matrix <- Gene_Resul$cov_matrix
# saveRDS(Gene_Cov_Matrix,glue::glue("{out_dir}/Mtb_gene_coverage.rds"))

#Get missing positions for each sample
missing_mat <- GetMissingMatrix(missing_files)
saveRDS(missing_mat,glue::glue("{out_dir}missing_mat.rds"))
