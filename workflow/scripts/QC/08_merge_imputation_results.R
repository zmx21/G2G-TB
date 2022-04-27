#Merge imputation results from AFGR and TB-DAR, based on better imputation accuracy per SNP.

library(glue)
library(dplyr)
library(data.table)

MergeImputationResults <- function(afgr_imputed_path,afgr_prefix,wgs_imputed_path,wgs_prefix,wgs_unimputed_info,out_dir,software_dir = '~/Software/',n_cores = 20){
  system(glue('mkdir -p {out_dir}'))
  system(glue('mkdir -p {out_dir}/tmp/'))
  
  #Read AFGR Imputed info 
  afgr_info <- data.table::fread(paste0(afgr_imputed_path,afgr_prefix,'.info.txt'),header = F)
  colnames(afgr_info) <- c('CHROM','POS','ID','REF','ALT','AF','INFO')
  afgr_nodup <- data.table::fread(paste0(afgr_imputed_path,afgr_prefix,'.nodup.vcf.gz.pos'),header = F)
  colnames(afgr_nodup) <- c('CHROM','POS','ID','REF','ALT')
  afgr_info <- dplyr::inner_join(afgr_nodup,afgr_info,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::select(CHROM,POS,ID=ID.x,REF,ALT,INFO,AF.AFGR=AF)
  remove(afgr_nodup)
  gc();
  
  #Extract WGS Imputed Info
  wgs_info <- rbind(data.table::fread(paste0(wgs_imputed_path,wgs_prefix,'.info.txt'),header = F))
  colnames(wgs_info) <- c('CHROM','POS','ID','REF','ALT','R2','AF')
  wgs_unimputed_info <- data.table::fread(wgs_unimputed_info)
  wgs_unimputed_info <- dplyr::select(wgs_unimputed_info,'CHROM'=V1,'POS'=V2,'ID'=V3,'REF'=V4,'ALT'=V5)
  wgs_info <- dplyr::inner_join(wgs_unimputed_info,wgs_info,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::select(CHROM,POS,ID.WGS=ID.x,ID.WGS.imp=ID.y,REF,ALT,R2,AF.WGS=AF)
  remove(wgs_unimputed_info)
  gc();
  
  #Merge WGS imputed and AFGR imputed
  wgs_afgr_imp_merged <- dplyr::full_join(afgr_info,wgs_info,by = c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::select(CHROM,POS,ID.AFGR=ID,REF,ALT,INFO,ID.WGS,ID.WGS.imp,R2,AF.WGS,AF.AFGR)
  afgr_snps <- which(wgs_afgr_imp_merged$INFO >= wgs_afgr_imp_merged$R2)
  afgr_snps <- union(afgr_snps,which(is.na(wgs_afgr_imp_merged$R2)))
  
  wgs_snps <- which(wgs_afgr_imp_merged$INFO < wgs_afgr_imp_merged$R2)
  wgs_snps <- union(wgs_snps,which(is.na(wgs_afgr_imp_merged$INFO)))
  
  write(wgs_afgr_imp_merged$ID.AFGR[afgr_snps],file=paste0(out_dir,'/tmp/afgr_imp_rsid.txt'))
  write(wgs_afgr_imp_merged$ID.WGS.imp[wgs_snps],file=paste0(out_dir,'/tmp/wgs_imp_rsid.txt'))
  write(wgs_afgr_imp_merged$ID.WGS[wgs_snps],file=paste0(out_dir,'/tmp/wgs_rsid.txt'))
  
  system(glue('{software_dir}plink2 --vcf {paste0(afgr_imputed_path,afgr_prefix,".nodup.vcf")} --extract {paste0(out_dir,"/tmp/afgr_imp_rsid.txt")} --export vcf bgz --out {paste0(out_dir,"/tmp/afgr_info_selected")}'))
  system(glue('{software_dir}plink2 --vcf {paste0(wgs_imputed_path,wgs_prefix,".vcf.gz")} --extract {paste0(out_dir,"/tmp/wgs_imp_rsid.txt")} --export vcf bgz --out {paste0(out_dir,"/tmp/wgs_imp_info_selected")}'))
  system(glue('{software_dir}bcftools index -t --threads {n_cores} {paste0(out_dir,"/tmp/wgs_imp_info_selected.vcf.gz")}'))
  system(glue('{software_dir}bcftools index -t --threads {n_cores} {paste0(out_dir,"/tmp/afgr_info_selected.vcf.gz")}'))

  #Merge Two reference panels
  system(glue('{software_dir}bcftools concat -a {paste0(out_dir,"/tmp/wgs_imp_info_selected.vcf.gz")} {paste0(out_dir,"/tmp/afgr_info_selected.vcf.gz")} -O z -o {paste0(out_dir,"/tmp/wgs.afgr.merged.vcf.gz")}'))
  system(glue('{software_dir}bcftools index -t --threads {n_cores} {paste0(out_dir,"/tmp/wgs.afgr.merged.vcf.gz")}'))
  system(glue('{software_dir}bcftools sort -T {out_dir} -O z -o {paste0(out_dir,"/tmp/wgs.afgr.merged.sorted.vcf.gz")} {paste0(out_dir,"/tmp/wgs.afgr.merged.vcf.gz")}'))
  system(glue('{software_dir}bcftools index -t --threads {n_cores} {paste0(out_dir,"/tmp/wgs.afgr.merged.sorted.vcf.gz")}'))
  
  #Write out merged file 
  wgs_afgr_imp_merged_afgr <- wgs_afgr_imp_merged[afgr_snps,] %>% dplyr::select(CHROM,POS,ID=ID.AFGR,REF,ALT,INFO,AF=AF.WGS,AF.AFGR)
  wgs_afgr_imp_merged_afgr$AF[is.na(wgs_afgr_imp_merged_afgr$AF)] <- wgs_afgr_imp_merged_afgr$AF.AFGR[is.na(wgs_afgr_imp_merged_afgr$AF)]
  wgs_afgr_imp_merged_afgr <- dplyr::select(wgs_afgr_imp_merged_afgr,-'AF.AFGR')
  wgs_afgr_imp_merged_wgs <- wgs_afgr_imp_merged[wgs_snps,] %>% dplyr::select(CHROM,POS,ID=ID.WGS,REF,ALT,INFO=R2,AF=AF.WGS)
  wgs_afgr_imp_merged_afgr$Source = 'AFGR'
  wgs_afgr_imp_merged_wgs$Source = 'WGS'
  wgs_afgr_imp_merged <- rbind(wgs_afgr_imp_merged_afgr,wgs_afgr_imp_merged_wgs)                                                     
  data.table::fwrite(wgs_afgr_imp_merged,file = paste0(out_dir,'WGS_AFGR_Merged_Info.txt'),sep = ' ')
}

afgr_imputed_path <- '../Imputed/AFGR/'
afgr_prefix <- 'TBDAR.AFGR.Imputed'
wgs_imputed_path <- '../Imputed/Tanz/'
wgs_prefix <- 'TBDAR.Tanz.Imputed'
out_dir <- '../Imputed/AFGR_Tanz/'
wgs_unimputed_info <- '../../../WGS/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz.info'

MergeImputationResults(afgr_imputed_path,afgr_prefix,wgs_imputed_path,wgs_prefix,wgs_unimputed_info,out_dir)