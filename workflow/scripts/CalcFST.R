library(dplyr)
Pop_Df <- data.table::fread('../../data/1000_Genomes/20130606_g1k.ped')
Pop_Info <- data.table::fread('../../data/1000_Genomes/20131219.populations.tsv') %>% dplyr::select(POP = `Population Code`,Super_POP = `Super Population`)
Pop_Df <- Pop_Df %>% dplyr::select(FID = `Individual ID`,IID = `Individual ID`,POP = Population) %>% dplyr::left_join(Pop_Info,by=c('POP' = 'POP')) %>% dplyr::filter(Super_POP !='AMR')

system('mkdir -p ../../results/FST/')
for(cur_pop in unique(Pop_Df$Super_POP)){
  data.table::fwrite(dplyr::filter(Pop_Df,Super_POP == cur_pop) %>% dplyr::select(IID),file = glue::glue('../../results/FST/{cur_pop}.txt'),col.names = F,row.names = F)

  # system(glue::glue("~/Software/bcftools view -S ../../results/FST/{cur_pop}.txt --force-samples ../../data/1000_Genomes/joined.1000genomes.vcf.gz -O v -o ../../results/FST/{cur_pop}.vcf" ))
  # system(glue::glue("~/Software/plink --vcf ../../results/FST/{cur_pop}.vcf --maf 0.05 --make-bed --out ../../results/FST/{cur_pop}.common.SNPs" ))
}

system(glue::glue("~/Software/bcftools view -r 18:71685316 ../../data/1000_Genomes/joined.1000genomes.vcf.gz -O z -o ../../results/FST/rs75769176.vcf.gz"))
system(glue::glue("~/Software/vcftools --gzvcf ../../results/FST/rs75769176.vcf.gz {paste(paste0('--weir-fst-pop ../../results/FST/',unique(Pop_Df$Super_POP),'.txt'),collapse = ' ')} --out ../../results/FST/rs75769176
"))

system(glue::glue("~/Software/bcftools view -r 21:43289997 ../../data/1000_Genomes/joined.1000genomes.vcf.gz -O z -o ../../results/FST/rs12151990.vcf.gz"))
system(glue::glue("~/Software/vcftools  --gzvcf ../../results/FST/rs12151990.vcf.gz {paste(paste0('--weir-fst-pop ../../results/FST/',unique(Pop_Df$Super_POP),'.txt'),collapse = ' ')} --out ../../results/FST/rs12151990
"))


system(glue::glue("~/Software/bcftools filter -e 'AFR_AF < 0.05 & EAS_AF < 0.05 & EUR_AF < 0.05 & SAS_AF < 0.05' ../../data/1000_Genomes/joined.1000genomes.vcf.gz --threads 20 -O z -o ../../results/FST/joined.1000genomes.commonSNPs.vcf.gz"))
system(glue::glue("~/Software/vcftools  --gzvcf ../../results/FST/joined.1000genomes.commonSNPs.vcf.gz {paste(paste0('--weir-fst-pop ../../results/FST/',unique(Pop_Df$Super_POP),'.txt'),collapse = ' ')} --out ../../results/FST/1KG_FST
"))

RASID_Chr18 <- data.table::fread('../../results/RASiD/RAiSD_Report.TBDAR.chr18',col.names = c('POS','mu')) 
Chr18_SNP_Pos <- which.min(abs(RASID_Chr18$POS - 71685316))
ggplot2::ggplot(data = RASID_Chr18 %>% dplyr::filter(POS > 6e07 & POS < 9e07)) + geom_point(aes(x=POS,y=mu),size = 0.1) + ggrepel::geom_label_repel(aes(x=POS,y=mu,label = 'rs75769176'),data = RASID_Chr18[Chr18_SNP_Pos,])  + geom_point(data = RASID_Chr18[Chr18_SNP_Pos,],aes(x=POS,y=mu),size = 1,color = 'red') + geom_hline(yintercept = quantile(RASID_Chr18$mu,0.9995)) + xlab('Chr18 POS')

RASID_Chr21 <- data.table::fread('../../results/RASiD/RAiSD_Report.TBDAR.chr21',col.names = c('POS','mu')) 
Chr21_SNP_Pos <- which.min(abs(RASID_Chr21$POS - 41869888))
ggplot2::ggplot(data = RASID_Chr18 %>% dplyr::filter(POS > 3e07 & POS < 6e07)) + geom_point(aes(x=POS,y=mu),size = 0.1) + ggrepel::geom_label_repel(aes(x=POS,y=mu,label = 'rs12151990'),data = RASID_Chr21[Chr21_SNP_Pos,])  + geom_point(data = RASID_Chr21[Chr21_SNP_Pos,],aes(x=POS,y=mu),size = 1,color = 'red') + geom_hline(yintercept = quantile(RASID_Chr21$mu,0.9995)) + xlab('Chr21 POS')

FST <- data.table::fread('../../results/FST/1KG_FST.weir.fst')
