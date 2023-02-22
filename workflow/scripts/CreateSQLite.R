library(DBI)
library(dplyr)
library(biomaRt)

CreateSQLLite <- function(){
  #Create DB if not already
  if(!file.exists("~/G2G_TB/results/SQL/TBDAR_G2G.sqlite")){
    db_con <- dbConnect(RSQLite::SQLite(), "~/G2G_TB/results/SQL/TBDAR_G2G.sqlite")
    dbDisconnect(db_con)
  }
  db_con <- dbConnect(RSQLite::SQLite(), "~/G2G_TB/results/SQL/TBDAR_G2G.sqlite")
  
  G2G_Dir <- '/home/zmxu/G2G_TB/results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/'
  G2G_Results <- dir(G2G_Dir)[sapply(dir(G2G_Dir),function(x) grepl(x=x,pattern = '.hybrid.gz'))]
  
  #Open each G2G Results file, append to SQL database
  for(i in 1:length(G2G_Results)){
    print(i)
    #Parse Mtb variant ID
    cur_Mtb_variant = paste0(strsplit(paste0(strsplit(G2G_Results[i],split = ':')[[1]][c(1,3)],collapse = ':'),'\\.')[[1]][1:2],collapse = '.')
    cur_DB <- data.table::fread(cmd = glue::glue("zcat {G2G_Dir}{G2G_Results[i]}"),select = c('ID','OR','P'))
    cur_DB$Mtb_Variant <- cur_Mtb_variant
    cur_DB <- cur_DB %>% dplyr::rename(Host_Variant = ID)
    #Append to SQL Table
    dbWriteTable(db_con,'TBDAR_G2G',cur_DB,append = T)
  }
  #Create Index on Host and Mtb Variant
  dbExecute(db_con, 'CREATE INDEX Host_ID ON TBDAR_G2G (Host_Variant,Mtb_Variant);')
  dbDisconnect(db_con)
  
}
FindTopAssociation <- function(Host_SNP_ID,Proxy_Type = NULL){
  print(Host_SNP_ID)
  db_con <- dbConnect(RSQLite::SQLite(), "~/G2G_TB/results/SQL/TBDAR_G2G.sqlite")
  Genotype_Path <- '../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/LINEAGE_ALL/TB_DAR_G2G'
  
  bim_file <- data.table::fread(glue::glue("{Genotype_Path}.bim"))
  
  if(Host_SNP_ID %in% bim_file$V2){
    #Query top G2G association with specified SNP
    query <- glue::glue("SELECT * FROM TBDAR_G2G WHERE \"Host_Variant\" = \"{Host_SNP_ID}\" ORDER BY \"P\" ASC LIMIT 1;")
    res <- dbGetQuery(db_con, query)
    dbDisconnect(db_con)
    
  }else if(is.null(Proxy_Type)){
    res <- data.frame(Host_Variant = Host_SNP_ID,OR = NA,P = NA,Mtb_Variant = NA)
    
    
  }else if(Proxy_Type == 'Region'){
    
    if(!file.exists(glue::glue("{Genotype_Path}.vcf.gz"))){
      system(glue::glue('~/Software/plink2 --bfile {Genotype_Path} --export vcf bgz --out {Genotype_Path}'))
      system(glue::glue('~/Software/bcftools index -t --threads 3 {Genotype_Path}.vcf.gz'))
    }
    snp_mart <- useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org",path="/biomart/martservice", dataset="hsapiens_snp")
    snp_pos <- getBM(attributes = c('chr_name', 'chrom_start', 'chrom_end'), 
                     filters = c('snp_filter'), 
                     values = Host_SNP_ID, 
                     mart = snp_mart)
    
    #Get SNPs within 10000 bp 
    bcf_query <- glue::glue("~/Software/bcftools view -r {snp_pos$chr_name}:{snp_pos$chrom_start - 5000}-{snp_pos$chrom_end + 5000} {Genotype_Path}.vcf.gz | ~/Software/bcftools query -f '%ID\n'")
    SNPs_to_incl <- system(bcf_query,intern = T)
    res <- data.frame(Host_Variant = character(),OR = numeric(),P = numeric(),Mtb_Variant = character())
    for(i in 1:length(SNPs_to_incl)){
      query <- glue::glue("SELECT * FROM TBDAR_G2G WHERE \"Host_Variant\" = \"{SNPs_to_incl[i]}\" ORDER BY \"P\" ASC LIMIT 1;")
      res <- rbind(res,dbGetQuery(db_con, query))
    }
    res <- res[which.min(res$P),]
  }
  
  res$Human_Query_SNP <- Host_SNP_ID
  res <- res %>% dplyr::rename(Human_Proxy_SNP = Host_Variant,Mtb_SNP = Mtb_Variant) %>% dplyr::relocate(Human_Query_SNP,Human_Proxy_SNP,Mtb_SNP,OR,P)
  return(res)
}
Susceptibility_SNPs <- c('rs2057178','rs4331426','rs10956514','rs4733781','rs12437118','rs6114027','rs73226617','rs34536443','rs17155120','rs4240897','rs4921437')
Susceptibility_Results <- lapply(Susceptibility_SNPs,function(x) FindTopAssociation(x,Proxy_Type = 'Region'))
data.table::fwrite(do.call(rbind,Susceptibility_Results),'../../results/Replication/Susceptibility.txt')

# Clade_SNPs <- c('rs17235409','rs114945555','rs3130660','rs529920','rs41472447')
# Clade_Results <- lapply(Clade_SNPs,function(x) FindTopAssociation(x,Proxy_Type = 'Region'))
# data.table::fwrite(do.call(rbind,Clade_Results),'../../results/Replication/Clade.txt')
# 
# Phelan_et_al_SNPs <- c('rs267951','rs74875032','rs60284130','rs142600697','rs1118438','rs558237','rs59441182','rs4563899')
# Phelan_et_al_Results <- lapply(Phelan_et_al_SNPs,function(x) FindTopAssociation(x,Proxy_Type = 'Region'))
# data.table::fwrite(do.call(rbind,Phelan_et_al_Results),'../../results/Replication/Phelan_et_al.txt')


