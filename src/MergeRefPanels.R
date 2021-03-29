library(dplyr)
library(ggrepel)

#Concanete sequenced and imputed samples, using common set of SNPs
ConcatSequencedImputed <- function(wgs_path,imputed_path,out_path){
  #Specify thresholds (both on single file and merged)
  af_diff_thresh <- 0.1 #Difference in AF between imputed and WGS dataset
  snp_missingness_thresh <- 0.1
  hwe_thresh <- 1e-6 #According to 1KG https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4750478/
  bcftools <- '/home/zmxu/Software/bcftools'
  KING <- '/home/zmxu/Software/king'
  PLINK <- '/home/zmxu/Software/plink'
  GCTA = '/home/zmxu/Software/gcta'
  n_cores <- 40
  
  #Create output dir
  system(glue::glue("mkdir -p {out_path}"))
  system(glue::glue("mkdir -p {out_path}tmp/"))
  system(glue::glue("mkdir -p {out_path}figs/"))
  
  #Convert Imputed VCF to PLINK
  # system(
  #   glue::glue(
  #     "{PLINK}2 --vcf {imputed_path} --threads {n_cores} --set-all-var-ids '@:#[b37]$r,$a' --keep-allele-order --make-bed --out {out_path}tmp/imputed.tmp1"
  #   )
  # )
  
  ##Filter WGS data for HWE and missingness
  # system(
  #   glue::glue(
  #     "{PLINK}2 --bfile {wgs_path} --threads {n_cores} --geno {snp_missingness_thresh} --keep-allele-order --set-all-var-ids '@:#[b37]$r,$a' --autosome --make-bed --out {out_path}tmp/wgs.tmp1"
  #   )
  # )
  # system(
  #   glue::glue(
  #     "{PLINK} --bfile {out_path}tmp/wgs.tmp1 --threads {n_cores} --keep-allele-order --list-duplicate-vars suppress-first --make-bed --out {out_path}tmp/wgs.tmp2"
  #   )
  # )
  # system(
  #   glue::glue(
  #     "{PLINK} --bfile {out_path}tmp/wgs.tmp1 --hwe {hwe_thresh} --threads {n_cores} --exclude {out_path}tmp/wgs.tmp2.dupvar --keep-allele-order --make-bed --out {out_path}tmp/wgs.tmp3"
  #   )
  # )
  
  ##Get consensus set of SNPs
  # imputed_snps <- data.table::fread(glue::glue("{out_path}tmp/imputed.tmp1.bim"))
  # wgs_snps <- data.table::fread(glue::glue("{out_path}tmp/wgs.tmp3.bim"))
  # 
  # consensus_exact_match <- dplyr::inner_join(imputed_snps,wgs_snps,c('V1'='V1','V4'='V4','V5'='V5','V6'='V6'))
  # write(consensus_exact_match$V2.x,file = file(glue::glue("{out_path}tmp/imputed.consensus")))
  # write(consensus_exact_match$V2.y,file = file(glue::glue("{out_path}tmp/wgs.consensus")))
  
  #Genotyped
  # system(
  #   glue::glue(
  #     "{PLINK} --bfile {out_path}tmp/imputed.tmp1 --extract {out_path}tmp/imputed.consensus --threads {n_cores} --keep-allele-order --make-bed --out {out_path}tmp/imputed.consensus"
  #   )
  # )
  # system(
  #   glue::glue(
  #     "{PLINK}2 --bfile {out_path}tmp/imputed.consensus --freq --out {out_path}tmp/imputed.consensus"
  #   )
  # )
  
  #WGS
  # system(
  #   glue::glue(
  #     "{PLINK} --bfile {out_path}tmp/wgs.tmp3 --extract {out_path}tmp/wgs.consensus --threads {n_cores} --keep-allele-order --make-bed --out {out_path}tmp/wgs.consensus"
  #   )
  # )
  # system(
  #   glue::glue(
  #     "{PLINK}2 --bfile {out_path}tmp/wgs.consensus --freq --out {out_path}tmp/wgs.consensus"
  #   )
  # )
  
  #Merge Imputed and WGS
  # system(
  #   glue::glue(
  #     "{PLINK} --bfile {out_path}tmp/imputed.consensus --bmerge {out_path}tmp/wgs.consensus --keep-allele-order --make-bed --out {out_path}/tmp/WGS.imputed.tmp"
  #   )
  # )
  
  # system(
  #   glue::glue(
  #     "{PLINK}2 --bfile {out_path}/tmp/WGS.imputed.tmp --freq --out {out_path}tmp/WGS.imputed.tmp"
  #   )
  # )
  
  # merged_freq <- data.table::fread(glue::glue("{out_path}tmp/WGS.imputed.tmp.afreq"))
  # merged_freq <- dplyr::filter(merged_freq,ALT_FREQS > 0.05 & ALT_FREQS < 0.95)
  # wgs_freq <- data.table::fread(glue::glue("{out_path}tmp/wgs.consensus.afreq"))
  # imputed_freq <- data.table::fread(glue::glue("{out_path}tmp/imputed.consensus.afreq"))
  # p <- ggplot2::ggplot(data.frame(dplyr::inner_join(imputed_freq %>% dplyr::filter(ID %in% merged_freq$ID),wgs_freq %>% dplyr::filter(ID %in% merged_freq$ID),c('ID'='ID')))) + aes(x = ALT_FREQS.x,y = ALT_FREQS.y) + geom_bin2d() +
  #   scale_fill_continuous(type = "viridis") +
  #   theme_bw() + ylab('Genotyped ALT AF') + xlab('WGS ALT AF')
  
  # freq_diff <- imputed_freq$ID[which(abs(imputed_freq$ALT_FREQS - wgs_freq$ALT_FREQS) > af_diff_thresh)]
  
  
  #MAF difference filtered file
  # write(freq_diff,file = file(glue::glue("{out_path}tmp/consensus.freq.diff")))

  # system(
  #   glue::glue(
  #     "{PLINK}2 --bfile {out_path}/tmp/WGS.imputed.tmp --exclude {out_path}tmp/consensus.freq.diff --export vcf bgz --keep-allele-order --out {out_path}WGS_AFGR.imputed"
  #   )
  # )
  
  # system(
  #   glue::glue(
  #     "{bcftools} index -t --threads {n_cores} {out_path}WGS_AFGR.imputed.vcf.gz"
  #   )
  # )
  
  
  #Write INFO score file 
  # Imputed_Info <- data.table::fread(glue::glue("{imputed_path}.info"))
  # Imputed_bim <- data.table::fread(glue::glue("{out_path}tmp/imputed.tmp1.bim"))
  # IDs_kept <- system(glue::glue("{bcftools} query -f '%ID\n' {out_path}WGS_AFGR.imputed.vcf.gz"),intern = T)
  # Info_kept <- Imputed_Info$V1[match(IDs_kept,Imputed_bim$V2)]
  # data.table::fwrite(data.frame(INFO = Info_kept),glue::glue("{out_path}WGS_AFGR.imputed.vcf.gz.info"),col.names = F,row.names = F,quote = F)
  
  #Calculate kingship
  
  # system(
  #   glue::glue(
  #     "{PLINK} --bfile {out_path}/tmp/WGS.imputed.tmp --exclude {out_path}tmp/consensus.freq.diff --make-bed --out {out_path}/tmp/WGS_AFGR.imputed"
  #   )
  # )
  # fam_file <- data.table::fread(glue::glue("{out_path}/tmp/WGS_AFGR.imputed.fam"))
  # fam_file$V1 <- fam_file$V2
  # data.table::fwrite(fam_file,quote = F,sep = ' ',col.names = F,row.names = F,file = glue::glue("{out_path}/tmp/WGS_AFGR.imputed.fam"))
  
  # system(
  #   glue::glue(
  #     "{KING} --related --degree 2 -b {out_path}/tmp/WGS_AFGR.imputed.bed --prefix {out_path}/WGS_AFGR.imputed"
  #   )
  # )
  
}
MergeReferencePanel <- function(primary_ref_panel,secondary_ref_panel,out_path,INFO = 0.8,n_cores = 20){
  #Create output dir
  system(glue::glue("mkdir -p {out_path}"))
  system(glue::glue("mkdir -p {out_path}tmp/"))
  bcftools <- '/home/zmxu/Software/bcftools'
  PLINK <- '/home/zmxu/Software/plink'

  primary_info <- data.table::fread(glue::glue("{gsub(x=primary_ref_panel,pattern = '.recode',replacement = '.vcf.gz.info')}"))
  primary_bim <- data.table::fread(glue::glue("{primary_ref_panel}.bim"))
  primary_ID <- primary_bim$V2
  primary_std_ID <- paste0(primary_bim$V1,':',primary_bim$V4,primary_bim$V5,',',primary_bim$V6)
  if('Rsq' %in% colnames(primary_info)){
    primary_info <- dplyr::rename(primary_info,V1 = Rsq)
    primary_info_to_incl <- as.vector(dplyr::filter(primary_info,V1 > INFO) %>% dplyr::select(SNP))
    primary_info_to_incl <- primary_std_ID[match(primary_info_to_incl$SNP,primary_ID)]
  }else{
    primary_info_to_incl <- primary_std_ID[which(primary_info$V1 > INFO)]
  }
  remove(primary_bim)

  secondary_info <- data.table::fread(glue::glue("{gsub(x=secondary_ref_panel,pattern = '.recode',replacement = '.vcf.gz.info')}"))
  secondary_bim <- data.table::fread(glue::glue("{secondary_ref_panel}.bim"))
  secondary_ID <- secondary_bim$V2
  secondary_std_ID <- paste0(secondary_bim$V1,':',secondary_bim$V4,secondary_bim$V5,',',secondary_bim$V6)

  if('Rsq' %in% colnames(secondary_info)){
    secondary_info <- dplyr::rename(secondary_info,V1 = Rsq)
    secondary_info_to_incl <- as.vector(dplyr::filter(secondary_info,V1 > INFO) %>% dplyr::select(SNP))
    secondary_info_to_incl <- secondary_std_ID[match(secondary_info_to_incl$SNP,secondary_ID)]
  }else{
    secondary_info_to_incl <- secondary_std_ID[which(secondary_info$V1 > INFO)]
  }
  remove(secondary_bim)

  primary_id_to_keep <- primary_ID[match(primary_info_to_incl,primary_std_ID)]
  write(primary_id_to_keep,file = glue::glue("{out_path}tmp/primary_id.txt"))
  secondary_id_to_keep <- secondary_ID[match(setdiff(secondary_info_to_incl,primary_info_to_incl),secondary_std_ID)]
  write(secondary_id_to_keep,file = glue::glue("{out_path}tmp/secondary_id.txt"))

  system(glue::glue("{PLINK}2 --bfile {primary_ref_panel} --threads {n_cores} --extract {out_path}tmp/primary_id.txt --keep-allele-order --export vcf-4.2 bgz --out {out_path}tmp/primary"))
  system(glue::glue("{PLINK}2 --bfile {secondary_ref_panel} --threads {n_cores} --extract {out_path}tmp/secondary_id.txt --keep-allele-order --export vcf-4.2 bgz --out {out_path}tmp/secondary"))
  system(glue::glue("{bcftools} index -t --threads 5 {out_path}tmp/primary.vcf.gz"))
  system(glue::glue("{bcftools} index -t --threads 5 {out_path}tmp/secondary.vcf.gz"))

  system(glue::glue("{bcftools} concat -a {out_path}tmp/primary.vcf.gz {out_path}tmp/secondary.vcf.gz | {bcftools} sort -T {out_path}tmp/ -O z -o {out_path}Tanz_AFGR.imputed.vcf.gz"))
  system(glue::glue("{bcftools} index -t --threads 5 {out_path}Tanz_AFGR.imputed.vcf.gz"))
  
  IDs <- system(glue::glue("{bcftools} query -f '%ID\n' {out_path}Tanz_AFGR.imputed.vcf.gz"),intern = T)
  all_info <- c(primary_info$V1[match(primary_id_to_keep,primary_ID)],secondary_info$V1[match(secondary_id_to_keep,secondary_ID)])
  info_scores <- all_info[match(IDs,c(primary_id_to_keep,secondary_id_to_keep))]
  data.table::fwrite(data.frame(INFO = info_scores),col.names = F,row.names = F,quote = F,glue::glue("{out_path}Tanz_AFGR.imputed.vcf.gz.info"))
  
}

MergeReferencePanel('/home/zmxu/G2G_TB/data/Genotyping/TB_DAR_AFGR_Imputed/AFGR.imputed.recode',
                    '/home/zmxu/G2G_TB/data/Genotyping/TB_DAR_Tanz_Imputed/Tanz.imputed.recode',
                    '/home/zmxu/G2G_TB/data/Genotyping/TB_DAR_AFGR_Tanz_Imputed/')
# ConcatSequencedImputed('/home/zmxu/G2G_TB/data/WGS/WGS_Host_Data/joined.hg19','/home/zmxu/G2G_TB/data/Genotyping/TB_DAR_AFGR_Imputed/AFGR.imputed.vcf.gz','/home/zmxu/G2G_TB/data/Genotyping/TB_DAR_WGS_AFGR_Imputed/')
