library(glue)
library(dplyr)
library(ggplot2)
library(ggrepel)
RunPCA1KG <- function(VCF_Path,Out_Path,maf_thresh = 0.01,hwe_thresh = 1e-6){
  #Set Path and variables
  n_cores <- 20
  PLINK = '/home/zmxu/Software/plink'
  GCTA = '/home/zmxu/Software/gcta'
  bcftools = '/home/zmxu/Software/bcftools'
  FILE_REF <- '~/G2G_TB/data/1000_Genomes/joined.1000genomes'
  excl_regions <- '~/G2G_TB/data/1000_Genomes/exclusion_regions_hg19.txt' #https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium
  super_pop <- data.table::fread('~/G2G_TB/data/1000_Genomes/20131219.populations.tsv')
  tbl_1KG <- data.table::fread('~/G2G_TB/data/1000_Genomes/20130606_g1k.ped') %>% dplyr::left_join(super_pop,by = c('Population' = 'Population Code'))
  
  ## TB_DAR ###
  #Remove Duplicate entries in VCF
  no_dup_vcf <- gsub(x=VCF_Path,pattern = '.vcf.gz',replacement = '.nodup.vcf.gz')
  system(
    glue::glue("{bcftools} norm -d snps {VCF_Path} | {bcftools} view -O z -o {no_dup_vcf}")
  )
  system(
    glue::glue("{bcftools} index -t {no_dup_vcf}")
  )
  
  bim_filt <- gsub(x=no_dup_vcf,pattern = '.vcf.gz',replacement = '.hwe')
  #Filter based on HWE and set SNP ID
  system(
    glue::glue(
      "{PLINK}2 --vcf {no_dup_vcf} --threads {n_cores} --hwe {hwe_thresh} --make-bed --keep-allele-order --set-all-var-ids @:#[b37]\\$r,\\$a  --out {bim_filt}"
    )
  )
  system(
    glue::glue(
      "{PLINK}2 --bfile {bim_filt} --threads {n_cores} --keep-allele-order --export vcf --const-fid --out {bim_filt}"
    )
  )
  
  # #Remove long range LD region
  system(
    glue::glue(
      "{PLINK}2 --bfile {bim_filt} --threads {n_cores} --exclude range {excl_regions} --make-bed --keep-allele-order --out {Out_Path}_tmp1 --autosome"
    )
  )

  ### 1KG ###
  #Solve SNP ID issue
  system(
    glue::glue(
      "{PLINK}2 --bfile {FILE_REF} --threads {n_cores} --make-bed --keep-allele-order --set-all-var-ids @:#[b37]\\$r,\\$a --out {Out_Path}_1KG --autosome"
    )
  )
  snps_1KG <- data.table::fread(glue::glue("{Out_Path}_1KG.bim"))
  snps_TB_DAR <- data.table::fread(glue::glue("{Out_Path}_tmp1.bim"))
  consensus_snps <- dplyr::inner_join(snps_TB_DAR,snps_1KG,by=c('V1'='V1','V4'='V4','V5'='V5','V6'='V6'))
  write(consensus_snps$V2.y,file = glue::glue("{Out_Path}_1KG_consensus_snps.txt"))
  system(
    glue::glue(
      "{PLINK}2 --bfile {Out_Path}_1KG --threads {n_cores} --keep-allele-order --make-bed --extract {Out_Path}_1KG_consensus_snps.txt --out {Out_Path}_1KG_consensus --autosome"
    )
  )

  #Generate TB-DAR with only consensus SNPs
  system(
    glue::glue(
      "{PLINK}2 --bfile {Out_Path}_tmp1 --threads {n_cores} --keep-allele-order --make-bed --extract {Out_Path}_1KG_consensus_snps.txt --out {Out_Path}_consensus --autosome"
    )
  )

  #Merge TB-DAR and 1KG
  system(
    glue::glue(
      "{PLINK} --bfile {Out_Path}_1KG_consensus --bmerge {Out_Path}_consensus --threads {n_cores} --keep-allele-order --make-bed --out {Out_Path}_1KG_TBDAR_consensus"
    )
  )

  #Prune (r2 < 0.2)
  system(
    glue::glue(
      "{PLINK} --bfile {Out_Path}_1KG_TBDAR_consensus --threads {n_cores} --keep-allele-order --indep-pairwise 200 100 0.2 --out {Out_Path}_1KG_TBDAR_consensus --autosome"
    )
  )

  system(
    glue::glue(
      "{PLINK} --bfile {Out_Path}_1KG_TBDAR_consensus --threads {n_cores} --keep-allele-order --extract {Out_Path}_1KG_TBDAR_consensus.prune.in --make-bed --out {Out_Path}_1KG_TBDAR_consensus_pruned --autosome"
    )
  )


  #Calculate PCA on the merged cohort
  n_pcs <- nrow(data.table::fread(glue::glue("{Out_Path}_1KG_TBDAR_consensus_pruned.fam")))
  system(
    glue::glue(
      "{PLINK} --bfile {Out_Path}_1KG_TBDAR_consensus_pruned --threads {n_cores} --pca {n_pcs} header var-wts --make-rel --out {Out_Path}_1KG_TBDAR_consensus_pruned"
    )
  )

  system(
    glue::glue(
      "{GCTA} --bfile {Out_Path}_1KG_TBDAR_consensus_pruned --make-grm --out {Out_Path}_1KG_TBDAR_consensus_pruned --thread-num {n_cores}"
    )
  )
  system(
    glue::glue(
      "{GCTA} --grm {Out_Path}_1KG_TBDAR_consensus_pruned --pca 20 --out {Out_Path}_1KG_TBDAR_consensus_pruned --thread-num {n_cores}"
    )
  )

  #Calculate PCA on the merged AFR cohort
  Non_AFR_Samples <- dplyr::filter(tbl_1KG,`Super Population` != 'AFR') %>% dplyr::select(`Individual ID`)
  consensus_fam_file <- data.table::fread(glue::glue("{Out_Path}_1KG_TBDAR_consensus.fam")) %>% dplyr::filter(V2 %in% Non_AFR_Samples$`Individual ID`) %>% dplyr::select(V1,V2)
  data.table::fwrite(consensus_fam_file,col.names = F,row.names = F,sep = ' ',file = glue::glue("{Out_Path}_1KG_Non_AFR.txt"))
  system(
    glue::glue(
      "{PLINK} --bfile {Out_Path}_1KG_TBDAR_consensus_pruned --remove {Out_Path}_1KG_Non_AFR.txt --keep-allele-order --make-bed --out {Out_Path}_1KG_AFR_TBDAR_consensus_pruned"
    )
  )

  system(
    glue::glue(
      "{GCTA} --bfile {Out_Path}_1KG_AFR_TBDAR_consensus_pruned --make-grm --out {Out_Path}_1KG_AFR_TBDAR_consensus_pruned --thread-num {n_cores}"
    )
  )
  system(
    glue::glue(
      "{GCTA} --grm {Out_Path}_1KG_AFR_TBDAR_consensus_pruned --pca 20 --out {Out_Path}_1KG_AFR_TBDAR_consensus_pruned --thread-num {n_cores}"
    )
  )

  #Cleanup
  system(glue::glue("rm {Out_Path}*tmp*"))

}
#Prepare 1KG sample table
VCF_Path <- '~/TB_GWAS/data/Genotyping/TB_DAR_Genotyping_WGS_Merged/WGS.genotyped.merged.vcf.gz'
Out_Path <- '~/TB_GWAS/PCA_WGS_Genotyped/TB_DAR_PCA'

# VCF_Path <- '~/TB_GWAS/data/Genotyping/TB_DAR_WGS_AFGR_Imputed/WGS_AFGR.imputed.nodup.vcf.gz'
# Out_Path <- '~/TB_GWAS/data/Genotyping/TB_DAR_WGS_AFGR_Imputed/PCA/TB_DAR_PCA'

# RunPCA1KG(VCF_Path,Out_Path)
super_pop <- data.table::fread('~/TB_array_design/1000_Genomes/20131219.populations.tsv')
tbl_1KG <- data.table::fread('~/TB_array_design/1000_Genomes/20130606_g1k.ped') %>% dplyr::left_join(super_pop,by = c('Population' = 'Population Code'))

#Load PCA results
pc_df <- data.table::fread(glue::glue("{Out_Path}_1KG_TBDAR_consensus_pruned.eigenvec"))
colnames(pc_df) <- c('FID','IID',paste0('PC',seq(1,ncol(pc_df) - 2)))
#Remove ID Prefix
pc_df$IID <- sapply(pc_df$IID,function(x) ifelse(grepl(pattern='_',x=x),strsplit(x=x,split = '_')[[1]][2],x))
#Plot with all Superpopulations
merged_df <- pc_df %>%
  dplyr::left_join(tbl_1KG,c('IID'='Individual ID'))
merged_df$`Data Set` <- as.factor(ifelse(is.na(merged_df$Population),'TB DAR','1000 Genomes'))
merged_df$`Super Population` <- as.factor(as.character(merged_df$`Super Population`))
merged_df$`Super Population` <-factor(merged_df$`Super Population` ,
                                      levels = rev(levels(merged_df$`Super Population`)))

eigen_val <- as.vector(data.table::fread(glue::glue("{Out_Path}_1KG_TBDAR_consensus_pruned.eigenval")))
var_explained <- eigen_val$V1 / sum(eigen_val$V1) * 100
merged_df$simple_ID <- sapply(merged_df$IID,function(x) ifelse(grepl(x=x,pattern = '-'),strsplit(x=x,split = '-')[[1]][1],x))
pc_plot_thousand_genome_PC1_PC2 <- ggplot(data = merged_df) +
  aes(x=PC1,y=PC2,color=`Super Population`,shape=`Data Set`) +
  geom_point() + scale_shape_manual(values = c(20,4)) +
  xlab(paste0('PC1 (',signif(var_explained[1],2),'%)'))  +
  ylab(paste0('PC2 (',signif(var_explained[2],2),'%)')) +
  geom_text_repel(data= dplyr::filter(merged_df,`Data Set` == 'TB DAR' & PC1 < 0),aes(label=simple_ID))
pc_plot_thousand_genome_PC3_PC4 <- ggplot(data = merged_df) +
  aes(x=PC3,y=PC4,color=`Super Population`,shape=`Data Set`) +
  geom_point() + scale_shape_manual(values = c(20,4)) +
  xlab(paste0('PC3 (',signif(var_explained[3],2),'%)'))  +
  ylab(paste0('PC4 (',signif(var_explained[4],2),'%)')) #+
pc_plot_thousand_genome <- ggpubr::ggarrange(plotlist = list(pc_plot_thousand_genome_PC1_PC2,pc_plot_thousand_genome_PC3_PC4),ncol = 2)

#AFR 1KG
# pc_df_AFR <- data.table::fread(glue::glue("{Out_Path}_1KG_AFR_TBDAR_consensus_pruned.eigenvec"))
# colnames(pc_df_AFR) <- c('FID','IID',paste0('PC',seq(1,ncol(pc_df_AFR) - 2)))
# #Remove ID Prefix
# pc_df_AFR$IID <- sapply(pc_df_AFR$IID,function(x) ifelse(grepl(pattern='_',x=x),strsplit(x=x,split = '_')[[1]][2],x))
# 
# merged_df_AFR <- pc_df_AFR %>% 
#   dplyr::left_join(tbl_1KG,c('IID'='Individual ID'))
# merged_df_AFR$`Data Set` <- as.factor(ifelse(is.na(merged_df_AFR$Population),'TB DAR','1000 Genomes'))
# merged_df_AFR$simple_ID <- sapply(merged_df_AFR$IID,function(x) ifelse(grepl(x=x,pattern = '-'),strsplit(x=x,split = '-')[[1]][1],x))
# 
# eigen_val_afr <- as.vector(data.table::fread(glue::glue("{Out_Path}_1KG_AFR_TBDAR_consensus_pruned.eigenval")))
# var_explained_afr <- eigen_val_afr / sum(eigen_val_afr) * 100
# 
# pc_plot_afr_pc1_pc2 <- ggplot(data = merged_df_AFR) + 
#   aes(x=PC1,y=PC2,color=`Population Description`,shape=`Data Set`) + 
#   geom_point() + scale_shape_manual(values = c(20,4)) +
#   xlab(paste0('PC1 (',signif(var_explained_afr[1],2),'%)'))  + 
#   ylab(paste0('PC2 (',signif(var_explained_afr[2],2),'%)')) +
#   geom_text_repel(data= dplyr::filter(merged_df_AFR,`Data Set` == 'TB DAR' & PC2 < -0.2),aes(label=simple_ID))
# pc_plot_afr_pc3_pc4 <- ggplot(data = merged_df_AFR) + 
#   aes(x=PC3,y=PC4,color=`Population Description`,shape=`Data Set`) + 
#   geom_point() + scale_shape_manual(values = c(20,4)) +
#   xlab(paste0('PC3 (',signif(var_explained_afr[3],2),'%)'))  + 
#   ylab(paste0('PC4 (',signif(var_explained_afr[4],2),'%)'))
# pc_plot_afr <- ggpubr::ggarrange(plotlist = list(pc_plot_afr_pc1_pc2,pc_plot_afr_pc3_pc4),ncol = 2)
# 
# pc_plots <- ggpubr::ggarrange(plotlist = list(pc_plot_thousand_genome_PC1_PC2,pc_plot_afr_pc1_pc2),ncol = 2)
