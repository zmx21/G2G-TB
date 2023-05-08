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
  ggrepel::geom_label_repel(mapping = aes(x=V2,label = 'Rv2348:p.Ile101Met', y = 0),data = dplyr::filter(lambda_values,V1 %in% c('Rv2348c_Rv2348c:2626678:p.Ile101Met')),nudge_y = 200) + 
  ggrepel::geom_label_repel(mapping = aes(x=V2,label = 'fixA:p.Thr67Met', y = 0),data = dplyr::filter(lambda_values,V1 %in% c('fixA_Rv3029c:3388671:p.Thr67Met')),nudge_y = 200) + ylim(0,200)
