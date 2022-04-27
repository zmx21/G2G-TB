library(networkD3)

system('~/Software/king -b ../../../Genotyping_WGS/TBDAR.WGS.Imputed.AnalysisReady.bed --related --degree 2 --cpus 50 --rplot --prefix ../../../Genotyping_WGS/TBDAR.WGS.Imputed.AnalysisReady')

king_dup_full <- data.table::fread('../../../Genotyping_WGS/TBDAR.WGS.Imputed.AnalysisReady.kin0') %>% 
  dplyr::filter(InfType == 'Dup/MZ') %>% dplyr::select(ID1,ID2)
king_dup_full$ID1_Simple <- sapply(king_dup_full$ID1,function(x) gsub(x=x,pattern = 'Batch1_',replacement = ''))
king_dup_full$ID1_Simple <- sapply(king_dup_full$ID1_Simple,function(x) gsub(x=x,pattern = 'Batch2_',replacement = ''))
king_dup_full$ID1_Simple <- sapply(king_dup_full$ID1_Simple,function(x) gsub(x=x,pattern = 'WGS_Fellay.',replacement = ''))

king_dup_full$ID2_Simple <- sapply(king_dup_full$ID2,function(x) gsub(x=x,pattern = 'Batch1_',replacement = ''))
king_dup_full$ID2_Simple <- sapply(king_dup_full$ID2_Simple,function(x) gsub(x=x,pattern = 'Batch2_',replacement = ''))
king_dup_full$ID2_Simple <- sapply(king_dup_full$ID2_Simple,function(x) gsub(x=x,pattern = 'WGS_Fellay.',replacement = ''))

genetic_duplicates <- dplyr::filter(king_dup_full,ID1_Simple != ID2_Simple)
genetic_duplicates <- c(genetic_duplicates$ID1,genetic_duplicates$ID2)
data.table::fwrite(data.frame(FID=genetic_duplicates,IID=genetic_duplicates),'../../../Genotyping_WGS/QC/Genetic_Duplicates.txt',sep = ' ',row.names = F,quote = F)

relatives <- data.table::fread('../../../Genotyping_WGS/TBDAR.WGS.Imputed.AnalysisReady.kin0') %>% 
  dplyr::filter(InfType != 'Dup/MZ') %>% dplyr::select(ID1,ID2,InfType)
relatives$ID1_Simple <- sapply(relatives$ID1,function(x) gsub(x=x,pattern = 'Batch1_',replacement = ''))
relatives$ID1_Simple <- sapply(relatives$ID1_Simple,function(x) gsub(x=x,pattern = 'Batch2_',replacement = ''))
relatives$ID1_Simple <- sapply(relatives$ID1_Simple,function(x) gsub(x=x,pattern = 'WGS_Fellay.',replacement = ''))

relatives$ID2_Simple <- sapply(relatives$ID2,function(x) gsub(x=x,pattern = 'Batch1_',replacement = ''))
relatives$ID2_Simple <- sapply(relatives$ID2_Simple,function(x) gsub(x=x,pattern = 'Batch2_',replacement = ''))
relatives$ID2_Simple <- sapply(relatives$ID2_Simple,function(x) gsub(x=x,pattern = 'WGS_Fellay.',replacement = ''))

relatives_pairs <- lapply(1:nrow(relatives),function(i) sort(c(relatives$ID1_Simple[i],relatives$ID2_Simple[i])))
relatives_pairs <- relatives_pairs[!duplicated(relatives_pairs)]
relatives_pairs_list <- data.frame(from = sapply(relatives_pairs,function(x) x[1]),to = sapply(relatives_pairs,function(x) x[2]))
#Manually check to see if there are trios 
p <- simpleNetwork(relatives_pairs_list, height="100px", width="100px")
trio <- c('82006','81995','80734')

#randomly exclude one of the pair, map back to PLINK_IDs. Include the parent from Trio 
relatives_to_excl <- c(relatives_pairs_list$from[!relatives_pairs_list$from %in% trio],'82006','80734')

#Map back to PLINK ID
relatives_to_excl_PLINK <- unique(unlist(sapply(relatives_to_excl,function(x) c(relatives$ID1,relatives$ID2)[grepl(x=c(relatives$ID1,relatives$ID2),pattern = x)])))
data.table::fwrite(data.frame(FID=relatives_to_excl_PLINK,IID=relatives_to_excl_PLINK),'../../../Genotyping_WGS/QC/Relatives.txt',sep = ' ')
