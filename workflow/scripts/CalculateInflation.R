library(GWAS.utils)
library(latex2exp)
# CalculateInflation <- function(GWAS_Path){
#   print(GWAS_Path)
# 
#   summary_stats <- data.table::fread(cmd = glue::glue('zcat {GWAS_Path}'),select = c('P'))
#   lambda <- GWAS.utils::genomic_inflation(P=summary_stats$P)
#   variant_name <- strsplit(GWAS_Path,split = "/")[[1]][length(strsplit(GWAS_Path,split = "/")[[1]])]
#   variant_name <- paste0(strsplit(variant_name,split = '\\.')[[1]][1:2],collapse = '.')
# 
#   data.table::fwrite(data.frame(variant_name=variant_name,lambda=lambda),gsub(x=GWAS_Path,pattern = '.gz',replacement = '.lambda'),sep = '\t',col.names = F)
#   data.table::fwrite(summary_stats,gsub(x=GWAS_Path,pattern = '.gz',replacement = '.pval'),col.names = F)
# 
# }
results_dir <- '../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/'
# GWAS_paths <- dir(results_dir)
# GWAS_paths <- GWAS_paths[grepl(GWAS_paths,pattern = '.gz')]
# lapply(GWAS_paths,function(x) CalculateInflation(glue::glue("{results_dir}{x}")))

lambda_paths <- dir(results_dir)
lambda_paths <- lambda_paths[grepl(lambda_paths,pattern = '.lambda')]

lambda_values <- do.call(rbind,lapply(lambda_paths,function(x) data.table::fread(glue::glue('{results_dir}{x}'))))

ggplot2::ggplot(data = lambda_values) + geom_histogram(aes(x=V2),bins = 40) + xlab(expression(paste("Genomic Inflation Factor (", lambda,')'))) + ylab('Count') +
  ggrepel::geom_text_repel(mapping = aes(x=V2),data = dplyr::filter(lambda_values,V1 %in% c('Rv2348c_Rv2348c:2626678:p.Ile101Met'))) 

# hist(lambda_values$V2,xlab = expression(lambda),main = expression(paste("Genomic Inflation Factor (", lambda,')')))

# G2G_Res <- readRDS('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/G2G_results.rds')
# G2G_Res_Simple <- lapply(G2G_Res$ALL,function(x) dplyr::select(x,'#CHROM',POS,P))
# 
# scale_factor <- 100
# host_chr_length <- getChromInfoFromUCSC("hg19") %>% dplyr::filter(chrom %in% c(paste0('chr',1:22),'chrX'))
# G2G_Df <- data.frame(Mtb_Pos=numeric(),Host_Pos=numeric(),P=numeric())
# for(i in 1:length(G2G_Res_Simple)){
#   cur_Mtb_pos <- as.numeric(strsplit(x=names(G2G_Res_Simple)[i],split = ':')[[1]][2])
#   for(j in 1:nrow(G2G_Res_Simple[[i]])){
#     if(nrow(G2G_Res_Simple[[i]]) == 0){
#       next
#     }
#     cur_chrom <- G2G_Res_Simple[[i]]$`#CHROM`[j]
#     if(paste0('chr',cur_chrom) == 'chr1'){
#       cum_length = 0
#     }else{
#       cum_length <- sum(host_chr_length$size[1:(which(host_chr_length$chrom == paste0('chr',cur_chrom)) - 1)]/ scale_factor)
#     }
#     cur_host_pos <- G2G_Res_Simple[[i]]$POS[j]/scale_factor + cum_length
#     cur_P <- G2G_Res_Simple[[i]]$P[j]
#     G2G_Df <- rbind(G2G_Df,data.frame(Mtb_Pos=cur_Mtb_pos,Host_Pos=cur_host_pos,P=cur_P))
#   }
# }
# 
# P_Thresh <- 1.02459e-10
# ggplot(data = G2G_Df) + geom_point(aes(x=Host_Pos,y=Mtb_Pos,color = -log10(P))) + 
#   scale_x_continuous(breaks = c(0,cumsum(host_chr_length$size / scale_factor)[1:22]),labels = host_chr_length$chrom,expand = c(0, 0),limits = c(0,max(G2G_Df$Host_Pos) + 200000)) + 
#   xlab('Human Genome Position') + ylab('M.tb Genome Position') + guides(shape='none') 
# 
# 
# 
