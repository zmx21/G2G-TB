library("dplyr")
library("ggplot2")
library("readr")
library("coloc")
library(susieR)
library("GenomicRanges")
library("seqminer")
library("AnnotationDbi")
library("org.Hs.eg.db")
library(LDlinkR)
library(cowplot)

GetGWASSummaryStats <- function(sum_stats_path,offset,maf_path = '../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/LINEAGE_ALL/TB_DAR_G2G.frq'){
  #Read MAF
  maf_file <- data.table::fread(maf_path) %>% dplyr::select(ID=SNP,MAF)
  
  #Import Summary stats
  sum_stats <- data.table::fread(cmd = glue::glue('zcat {sum_stats_path}')) %>% dplyr::select(CHR='#CHROM',POS,rsid=ID,REF,ALT,OR,SE='LOG(OR)_SE',P,OBS_CT)
  sum_stats_hits <- dplyr::filter(sum_stats,CHR == sum_stats$CHR[which.min(sum_stats$P)]) %>% dplyr::filter(POS < sum_stats$POS[which.min(sum_stats$P)] + offset & POS > sum_stats$POS[which.min(sum_stats$P)] - offset) 
  sum_stats_hits$CHR <- paste0('chr',sum_stats_hits$CHR)
  #Import chain file
  chain = rtracklayer::import.chain("../../data/hg19ToHg38.over.chain")
  
  #Convert to granges
  sum_stats_hits_granges <- makeGRangesFromDataFrame(sum_stats_hits %>% dplyr::left_join(maf_file,by=c('rsid'='ID')) ,keep.extra.columns=T,ignore.strand=T,seqnames.field	
= 'CHR',start.field	= 'POS',end.field	= 'POS')
  
  #Lift over summary statistics
  gwas_stats_hg38 = rtracklayer::liftOver(sum_stats_hits_granges, chain) %>% 
    unlist() %>% 
    dplyr::as_tibble() %>%
    dplyr::select(CHR = seqnames, POS = start,rsid,REF,ALT,OR,SE,P,MAF,OBS_CT)%>% 
    dplyr::mutate(beta = log(OR)) %>%
    dplyr::select(-OR) %>%
    dplyr::mutate(id = paste(CHR,POS,sep = ':'))  %>%   
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1)%>% 
    dplyr::mutate(id = paste0(gsub(CHR,pattern = 'chr',replacement = ''),':',POS)) 
    #dplyr::mutate(id = paste0(gsub(CHR,pattern = 'chr',replacement = ''),':',POS,'_',REF,'_',ALT)) 
  
  
}

GetEQTLSummaryStats <- function(GWAS_Sum_Stats,target_study,target_qtl_group,ensid){

  import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE){
    
    #Fetch summary statistics with seqminer
    fetch_table = seqminer::tabix.read.table(tabixFile = ftp_path, tabixRange = region, stringsAsFactors = FALSE) %>%
      dplyr::as_tibble()
    colnames(fetch_table) = column_names
    
    #Remove rsid duplicates and multi-allelic variant
    summary_stats = dplyr::filter(fetch_table, gene_id == selected_gene_id) %>% #rsid duplicates
      dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
      dplyr::group_by(id) %>% 
      dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
      dplyr::filter(row_count == 1) %>% 
      dplyr::mutate(id = paste0(chromosome,':',position))
    #dplyr::mutate(id = paste0(chromosome,':',position,'_',ref,'_',alt))
    return(summary_stats)
  }
  
  tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

  eqtl_region = paste0(gsub(unique(GWAS_Sum_Stats$CHR),pattern = 'chr',replacement = ''),':',min(GWAS_Sum_Stats$POS),'-',max(GWAS_Sum_Stats$POS))
  eqtl_df = dplyr::filter(tabix_paths, qtl_group == target_qtl_group & study == target_study) %>% 
    dplyr::filter(quant_method == "ge" | quant_method == "microarray")
  
  #Extract column names from first file
  column_names = colnames(readr::read_tsv(gsub(eqtl_df$ftp_path,pattern = 'ftp://',replacement = '~/G2G_TB/data/eQTL/'), n_max = 1))
  
  #Import summary statistics
  summary_stats = import_eQTLCatalogue(gsub(eqtl_df$ftp_path,pattern = 'ftp://',replacement = '~/G2G_TB/data/eQTL/'), eqtl_region, selected_gene_id = ensid, column_names)
  return(summary_stats)
}

GetLDMatrix <- function(snps,POP){
  LD_Mat <- matrix(NA,nrow = length(snps),ncol = length(snps))
  rownames(LD_Mat) <- snps
  colnames(LD_Mat) <- snps
  if(ncol(LD_Mat) > 500){
    row_bound <- c(seq(from = 1,to = nrow(LD_Mat),by = 499),nrow(LD_Mat))
    col_bound <- c(seq(from = 1,to = ncol(LD_Mat),by = 499),ncol(LD_Mat))
  }else{
    row_bound <- c(1,nrow(LD_Mat))
    col_bound <- c(1,ncol(LD_Mat))
  }

  for(i in 2:length(row_bound)){
    rows <- rownames(LD_Mat)[row_bound[i-1]:row_bound[i]]
    for(j in 2:length(col_bound)){
      cols <- colnames(LD_Mat)[col_bound[j-1]:col_bound[j]]
      cur_LD <- LDmatrix(snps = unique(c(rows,cols)), 
                         pop = POP, 
                         r2d = "r2", 
                         token = 'b9f8399d9bab',
                         genome_build = "grch38")
      cur_LD_rownames <- cur_LD[,1,drop=T]
      cur_LD <- as.matrix(cur_LD[,-1])
      rownames(cur_LD) <- cur_LD_rownames
      LD_Mat[intersect(rows,rownames(cur_LD)),intersect(cols,colnames(cur_LD))] <- cur_LD[intersect(rows,rownames(cur_LD)),intersect(cols,colnames(cur_LD))]
    }
  }
  return(LD_Mat)
}

RunColoc <- function(eqtl_sumstats,gwas_sumstats,out_path,eQTL_POP = 'EUR',method = 'Coloc',Pathogen_MAC){
  if(method == 'Susie' | method == 'Coloc_Cond' | method == 'Coloc_Mask'){
    system(glue::glue("mkdir -p {out_path}/tmp/"))
    EQTL_SNP_List <- eqtl_sumstats$variant
    write(EQTL_SNP_List,file = glue::glue("{out_path}/tmp/EQTL_SNPs"))
    Sample_Manifest_1KG <- dplyr::left_join(data.table::fread('../../data/1000_Genomes/20130606_g1k.ped'),
                                            data.table::fread('../../data/1000_Genomes/20131219.populations.tsv') %>% 
                                              dplyr::select(Population = `Population Code`,SuperPOP = `Super Population`),by =c('Population'='Population'))
    EQTL_Samples <- dplyr::filter(Sample_Manifest_1KG,SuperPOP %in% eQTL_POP)$`Individual ID`
    data.table::fwrite(data.frame(FID=0,IID=EQTL_Samples),file = glue::glue("{out_path}/tmp/EQTL_Samples"),row.names = F,col.names = F,sep = '\t')
    
    bfile_path_1KG <- glue::glue('../../data/1000_Genomes/GRCh38/chr{unique(eqtl_sumstats$chromosome)}.recalibrated_variants')
    system(glue::glue("/home/zmxu/Software/plink --bfile {bfile_path_1KG} --keep-allele-order --extract {out_path}/tmp/EQTL_SNPs --keep {out_path}/tmp/EQTL_Samples --make-bed --out {out_path}/tmp/EQTL_LD_BFILE"))
    system(glue::glue("/home/zmxu/Software/plink --bfile {out_path}/tmp/EQTL_LD_BFILE --keep-allele-order --r square --out {out_path}EQTL_LD"))
    
    
    eQTL_LD <- as.matrix(read.table(glue::glue("{out_path}EQTL_LD.ld")))
    rownames(eQTL_LD) <- data.table::fread(glue::glue("{out_path}/tmp/EQTL_LD_BFILE.bim"))$V2
    colnames(eQTL_LD) <- data.table::fread(glue::glue("{out_path}/tmp/EQTL_LD_BFILE.bim"))$V2
    
    eQTL_LD <- eQTL_LD[,!apply(eQTL_LD,2,function(x) all(is.nan(x)))]
    eQTL_LD <- eQTL_LD[!apply(eQTL_LD,1,function(x) all(is.nan(x))),]
    
    eQTL_LD[is.nan(eQTL_LD)] <- 0

    eqtl_sumstats <- eqtl_sumstats %>% dplyr::filter(variant %in% rownames(eQTL_LD))
    
    eQTL_dataset = list(beta = eqtl_sumstats$beta,
                        varbeta = eqtl_sumstats$se^2,
                        N = (eqtl_sumstats$an)[1]/2, # Samples size is allele number (AN) divided by 2
                        MAF = eqtl_sumstats$maf,
                        type = "quant",
                        snp = eqtl_sumstats$variant,
                        position = eqtl_sumstats$position,
                        LD = eQTL_LD)
    
    bfile_path <- '../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/LINEAGE_ALL/TB_DAR_G2G'
    gwas_sumstats$variant <- paste0(gwas_sumstats$CHR,'_',gwas_sumstats$POS,'_',gwas_sumstats$REF,'_',gwas_sumstats$ALT)

    write(gwas_sumstats$rsid,file = glue::glue("{out_path}/tmp/GWAS_SNPs"))
    system(glue::glue("/home/zmxu/Software/plink --bfile {bfile_path} --extract {out_path}/tmp/GWAS_SNPs --make-bed --keep-allele-order --out {out_path}/tmp/GWAS_LD_BFILE"))
    system(glue::glue("/home/zmxu/Software/plink --bfile {out_path}/tmp/GWAS_LD_BFILE --r square --keep-allele-order --out {out_path}GWAS_LD"))

    GWAS_LD <- as.matrix(read.table(glue::glue("{out_path}GWAS_LD.ld")))
    rownames(GWAS_LD) <- data.table::fread(glue::glue("{out_path}/tmp/GWAS_LD_BFILE.bim"))$V2
    colnames(GWAS_LD) <- data.table::fread(glue::glue("{out_path}/tmp/GWAS_LD_BFILE.bim"))$V2
    gwas_sumstats <- gwas_sumstats %>% dplyr::filter(rsid %in% rownames(GWAS_LD))
    rownames(GWAS_LD) <- gwas_sumstats$variant[match(rownames(GWAS_LD),gwas_sumstats$rsid)]
    colnames(GWAS_LD) <- gwas_sumstats$variant[match(colnames(GWAS_LD),gwas_sumstats$rsid)]
    
    
    # gwas_sumstats$variant <- paste0(gwas_sumstats$CHR,'_',gwas_sumstats$POS,'_',gwas_sumstats$REF,'_',gwas_sumstats$ALT)
    # write(gwas_sumstats$variant,file = glue::glue("{out_path}/tmp/GWAS_SNPs"))
    # Sample_Manifest_1KG <- dplyr::left_join(data.table::fread('../../data/1000_Genomes/20130606_g1k.ped'),
    #                                         data.table::fread('../../data/1000_Genomes/20131219.populations.tsv') %>%
    #                                           dplyr::select(Population = `Population Code`,SuperPOP = `Super Population`),by =c('Population'='Population'))
    # GWAS_Samples <- dplyr::filter(Sample_Manifest_1KG,Population == 'LWK')$`Individual ID`
    # data.table::fwrite(data.frame(FID=0,IID=GWAS_Samples),file = glue::glue("{out_path}/tmp/GWAS_Samples"),row.names = F,col.names = F,sep = '\t')
    # 
    # bfile_path_1KG <- glue::glue('../../data/1000_Genomes/GRCh38/chr{unique(eqtl_sumstats$chromosome)}.recalibrated_variants')
    # system(glue::glue("/home/zmxu/Software/plink --bfile {bfile_path_1KG} --keep-allele-order --extract {out_path}/tmp/GWAS_SNPs --keep {out_path}/tmp/GWAS_Samples --make-bed --out {out_path}/tmp/GWAS_LD_BFILE"))
    # system(glue::glue("/home/zmxu/Software/plink --bfile {out_path}/tmp/GWAS_LD_BFILE --keep-allele-order --r square --out {out_path}GWAS_LD"))
    # 
    # 
    # GWAS_LD <- as.matrix(read.table(glue::glue("{out_path}GWAS_LD.ld")))
    # rownames(GWAS_LD) <- data.table::fread(glue::glue("{out_path}/tmp/GWAS_LD_BFILE.bim"))$V2
    # colnames(GWAS_LD) <- data.table::fread(glue::glue("{out_path}/tmp/GWAS_LD_BFILE.bim"))$V2
    # 
    # GWAS_LD <- GWAS_LD[,!apply(GWAS_LD,2,function(x) all(is.nan(x)))]
    # GWAS_LD <- GWAS_LD[!apply(GWAS_LD,1,function(x) all(is.nan(x))),]
    # 
    # GWAS_LD[is.nan(GWAS_LD)] <- 0
    # 
    # gwas_sumstats <- gwas_sumstats %>% dplyr::filter(variant %in% rownames(GWAS_LD))
    
    gwas_dataset = list(beta = gwas_sumstats$beta,
                        varbeta = gwas_sumstats$SE^2,
                        type = "cc",
                        snp = gwas_sumstats$variant,
                        MAF = gwas_sumstats$MAF,
                        N = gwas_sumstats$OBS_CT,
                        position = gwas_sumstats$POS,
                        s = Pathogen_MAC/gwas_sumstats$OBS_CT,
                        LD = GWAS_LD)
    
    gwas_df <- gwas_sumstats %>% dplyr::select(rsid=variant,chr=CHR,pos=POS,pval=P)
    gwas_ld_df <- as.data.frame.table(GWAS_LD^2)
    colnames(gwas_ld_df) <- c('SNP_A','SNP_B','R2')

    eqtl_df <- eqtl_sumstats %>% dplyr::select(rsid=variant,chr=chromosome,pos=position,pval=pvalue)
    eQTL_ld_df <- as.data.frame.table(eQTL_LD^2)
    colnames(eQTL_ld_df) <- c('SNP_A','SNP_B','R2')
    source('./locuscomparer.R')
    p_locus <- locuscompare(in_fn1=eqtl_df,in_fn2=gwas_df,ld_fn1=eQTL_ld_df,ld_fn2=gwas_ld_df,combine=T,legend_position = 'topright')
    
    if(method == 'Susie'){
      eQTL_Susie <- runsusie(eQTL_dataset,coverage=0.95)
      GWAS_Susie <- runsusie(gwas_dataset,coverage=0.95,estimate_prior_variance=FALSE)
      
      susie.res <- NULL
      if(!is.null(eQTL_Susie)){
        if(requireNamespace("susieR",quietly=TRUE)) {
          susie.res=coloc.susie(eQTL_Susie,GWAS_Susie)
        }
      }
      saveRDS(list(result=susie.res,p_locus=p_locus,
                   gwas_dataset=gwas_dataset,
                   eQTL_dataset=eQTL_dataset,
                   GWAS_Susie=GWAS_Susie,
                   eQTL_Susie=eQTL_Susie),file = glue::glue("{out_path}Susie.rds"))
      if(is.null(susie.res)){
        PP <- rep(NA,5)
      }else{
        PP <- c(max(susie.res$summary$PP.H0.abf),max(susie.res$summary$PP.H1.abf),max(susie.res$summary$PP.H2.abf),max(susie.res$summary$PP.H3.abf),max(susie.res$summary$PP.H4.abf))
      }
      names(PP) <- c('PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')
      return(PP)
    }else if (method == 'Coloc_Cond'){
      gwas_dataset$method <- 'single'
      eQTL_dataset$method <- 'cond'
      coloc_cond <- coloc.signals(gwas_dataset,eQTL_dataset,p12 = 1e-6)
      PP <- coloc_cond$summary[,c('PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')]
      saveRDS(list(result=coloc_cond,p_locus=p_locus,
                   gwas_dataset=gwas_dataset,
                   eQTL_dataset=eQTL_dataset),file = glue::glue("{out_path}Coloc_Cond.rds"))
      
      
    }else if(method == 'Coloc_Mask'){
      gwas_dataset$method <- 'single'
      eQTL_dataset$method <- 'mask'
      coloc_mask <- coloc.signals(gwas_dataset,eQTL_dataset,p12 = 1e-6)
      PP <- coloc_mask$summary[,c('PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')]
      saveRDS(list(result=coloc_mask,p_locus=p_locus,
                   gwas_dataset=gwas_dataset,
                   eQTL_dataset=eQTL_dataset),file = glue::glue("{out_path}Coloc_Mask.rds"))
      return(PP)
    }

  }else if(method == 'Coloc'){
    gwas_sumstats$variant <- paste0(gwas_sumstats$CHR,'_',gwas_sumstats$POS,'_',gwas_sumstats$REF,'_',gwas_sumstats$ALT)
    
    gwas_dataset = list(beta = gwas_sumstats$beta,
                        varbeta = gwas_sumstats$SE^2,
                        type = "cc",
                        snp = gwas_sumstats$variant,
                        MAF = gwas_sumstats$MAF,
                        N = gwas_sumstats$OBS_CT,
                        position = gwas_sumstats$POS)
    
    eQTL_dataset = list(beta = eqtl_sumstats$beta,
                        varbeta = eqtl_sumstats$se^2,
                        N = (eqtl_sumstats$an)[1]/2, # Samples size is allele number (AN) divided by 2
                        MAF = eqtl_sumstats$maf,
                        type = "quant",
                        snp = eqtl_sumstats$variant,
                        position = eqtl_sumstats$position)
    coloc.res <- tryCatch(coloc.abf(gwas_dataset,eQTL_dataset),error = function(e) return(NULL))
    
    if(!is.null(coloc.res)){
      PP <- coloc.res$summary[-1]
    }else{
      PP <- rep(NA,5)
      names(PP) <- c('PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')
    }

    saveRDS(list(result=coloc.res,
                 gwas_dataset=gwas_dataset,
                 eQTL_dataset=eQTL_dataset),file = glue::glue("{out_path}Coloc.rds"))

    return(PP)
  }
}

RunAll <- function(GWAS_Path,Target_Gene,method = 'Coloc',Pathogen_MAC){
  offset <- 300000
  GWAS_stats <- GetGWASSummaryStats(GWAS_Path,offset)
  
  tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    dplyr::filter(quant_method == 'ge' | quant_method == 'microarray') %>% dplyr::filter(study == 'Quach_2016' & qtl_group == 'monocyte_naive')
  
  #Get eGene ENSEMBLE ID
  ensid = mapIds(org.Hs.eg.db,
                 keys=Target_Gene, 
                 column="ENSEMBL",
                 keytype="SYMBOL",
                 multiVals="first")
  
  PP <- lapply(1:nrow(tabix_paths),function(i){
    EQTL_stats <- GetEQTLSummaryStats(GWAS_stats,target_study = tabix_paths$study[i],target_qtl_group = tabix_paths$qtl_group[i],ensid = ensid) %>% 
      dplyr::filter(rsid != "NA\r")
    
    if(nrow(EQTL_stats) == 0){
      PP <- rep(NA,5)
      names(PP) <- c('PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')
      return(PP)
    }
    EQTL_stats$rsid <- gsub(EQTL_stats$rsid,pattern = '\r',replacement = '')
    
    out_path <- glue::glue('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/{method}/{Target_Gene}/{tabix_paths$study[i]}/{tabix_paths$qtl_group[i]}/')
    system(glue::glue('mkdir -p {out_path}'))
    return(RunColoc(EQTL_stats,GWAS_stats,out_path,method = method,Pathogen_MAC=Pathogen_MAC))
  })
  
  names(PP) <- paste0(tabix_paths$study,':',tabix_paths$qtl_group)
  saveRDS(PP,file = glue::glue('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/{method}/{Target_Gene}/PP.rds'))
  
}
RunAll('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/Rv2348c_Rv2348c:2626678:p.Ile101Met.Rv2348c_Rv2348c:2626678:p.Ile101Met.glm.logistic.hybrid.gz','PRDM15',method = 'Coloc_Mask',Pathogen_MAC = 18)
# RunAll('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/Rv2348c_Rv2348c:2626678:p.Ile101Met.Rv2348c_Rv2348c:2626678:p.Ile101Met.glm.logistic.hybrid.gz','PRDM15',method = 'Coloc_Cond',Pathogen_MAC = 18)

# RunAll('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/Rv2348c_Rv2348c:2626678:p.Ile101Met.Rv2348c_Rv2348c:2626678:p.Ile101Met.glm.logistic.hybrid.gz','PRDM15',method = 'Coloc',Pathogen_MAC = 18)
# RunAll('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/fixA_Rv3029c:3388671:p.Thr67Met.fixA_Rv3029c:3388671:p.Thr67Met.glm.logistic.hybrid.gz','FBXO15',method = 'Coloc')
# RunAll('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/fixA_Rv3029c:3388671:p.Thr67Met.fixA_Rv3029c:3388671:p.Thr67Met.glm.logistic.hybrid.gz','TIMM21',method = 'Coloc')
# RunAll('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/Rv2348c_Rv2348c:2626678:p.Ile101Met.Rv2348c_Rv2348c:2626678:p.Ile101Met.glm.logistic.hybrid.gz','C2CD2',method = 'Coloc')
# RunAll('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/Rv2348c_Rv2348c:2626678:p.Ile101Met.Rv2348c_Rv2348c:2626678:p.Ile101Met.glm.logistic.hybrid.gz','RIPK4',method = 'Coloc')

# RunAll('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/Rv2348c_Rv2348c:2626678:p.Ile101Met.Rv2348c_Rv2348c:2626678:p.Ile101Met.glm.logistic.hybrid.gz','PRDM15',method = 'Susie',Pathogen_MAC = 18)
# RunAll('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/fixA_Rv3029c:3388671:p.Thr67Met.fixA_Rv3029c:3388671:p.Thr67Met.glm.logistic.hybrid.gz','FBXO15',method = 'Susie')
# RunAll('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/fixA_Rv3029c:3388671:p.Thr67Met.fixA_Rv3029c:3388671:p.Thr67Met.glm.logistic.hybrid.gz','TIMM21',method = 'Susie')
# RunAll('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/Rv2348c_Rv2348c:2626678:p.Ile101Met.Rv2348c_Rv2348c:2626678:p.Ile101Met.glm.logistic.hybrid.gz','C2CD2',method = 'Susie',Pathogen_MAC = 18)
# RunAll('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/Rv2348c_Rv2348c:2626678:p.Ile101Met.Rv2348c_Rv2348c:2626678:p.Ile101Met.glm.logistic.hybrid.gz','RIPK4',method = 'Susie',Pathogen_MAC = 18)

# PRDM15_Coloc <- readRDS('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/Coloc/PRDM15/PP.rds')
# PRDM15_Coloc_H4 <- sapply(PRDM15_Coloc,function(x) x["PP.H4.abf"])
# PRDM15_Coloc_H4_Df <- data.frame(Tissue = names(PRDM15_Coloc),PP_H4 = PRDM15_Coloc_H4,Category = sapply(names(PRDM15_Coloc_H4),function(x) strsplit(strsplit(x=x,split = ':')[[1]][2],split = '_')[[1]][1]))
# PRDM15_Coloc_H4_Df$Category[!PRDM15_Coloc_H4_Df$Category %in% names(table(PRDM15_Coloc_H4_Df$Category))[table(PRDM15_Coloc_H4_Df$Category) > 2]] <- 'Other'
# 
# ggplot2::ggplot(data = PRDM15_Coloc_H4_Df %>% dplyr::filter(PP_H4 > 0.05),aes(x=reorder(Tissue,-PP_H4),y=PP_H4)) + geom_bar(stat = 'identity') + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(~Category,scales = 'free_x',nrow=2) + ylim(0,1) + ylab('Posterior - Shared Causal Variant') + xlab('Tissue') + ggtitle('Coloc - PRDM15~rs1215990')
# 
# C2CD2_Coloc <- readRDS('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/Coloc/C2CD2/PP.rds')
# C2CD2_Coloc_H4 <- sapply(C2CD2_Coloc,function(x) x["PP.H4.abf"])
# C2CD2_Coloc_H4_Df <- data.frame(Tissue = names(C2CD2_Coloc),PP_H4 = C2CD2_Coloc_H4,Category = sapply(names(C2CD2_Coloc),function(x) strsplit(strsplit(x=x,split = ':')[[1]][2],split = '_')[[1]][1]))
# C2CD2_Coloc_H4_Df$Category[!C2CD2_Coloc_H4_Df$Category %in% names(table(C2CD2_Coloc_H4_Df$Category))[table(C2CD2_Coloc_H4_Df$Category) > 2]] <- 'Other'
# 
# ggplot2::ggplot(data = C2CD2_Coloc_H4_Df %>% dplyr::filter(PP_H4 > 0.05),aes(x=reorder(Tissue,-PP_H4),y=PP_H4)) + geom_bar(stat = 'identity') + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(~Category,scales = 'free_x',nrow=2) + ylim(0,1) + ylab('Posterior - Shared Causal Variant') + xlab('Tissue') + ggtitle('Coloc - C2CD2~rs1215990')
# 
# FBXO15_Coloc <- readRDS('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/Coloc/FBXO15/PP.rds')
# FBXO15_Coloc_H4 <- sapply(FBXO15_Coloc,function(x) x["PP.H4.abf"])
# FBXO15_Coloc_H4_Df <- data.frame(Tissue = names(FBXO15_Coloc),PP_H4 = FBXO15_Coloc_H4,Category = sapply(names(FBXO15_Coloc),function(x) strsplit(strsplit(x=x,split = ':')[[1]][2],split = '_')[[1]][1]))
# FBXO15_Coloc_H4_Df$Category[!FBXO15_Coloc_H4_Df$Category %in% names(table(FBXO15_Coloc_H4_Df$Category))[table(FBXO15_Coloc_H4_Df$Category) > 2]] <- 'Other'
# 
# ggplot2::ggplot(data = FBXO15_Coloc_H4_Df %>% dplyr::filter(PP_H4 > 0.05),aes(x=reorder(Tissue,-PP_H4),y=PP_H4)) + geom_bar(stat = 'identity') + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + facet_wrap(~Category,scales = 'free_x',nrow=1) + ylim(0,1) + ylab('Posterior - Shared Causal Variant') + xlab('Tissue') + ggtitle('Coloc - FBXO15~rs1215990')
# 
# TIMM21_Coloc <- readRDS('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/Coloc/TIMM21/PP.rds')
# TIMM21_Coloc_H4 <- sapply(TIMM21_Coloc,function(x) x["PP.H4.abf"])
# TIMM21_Coloc_H4_Df <- data.frame(Tissue = names(TIMM21_Coloc),PP_H4 = TIMM21_Coloc_H4,Category = sapply(names(TIMM21_Coloc),function(x) strsplit(strsplit(x=x,split = ':')[[1]][2],split = '_')[[1]][1]))
# TIMM21_Coloc_H4_Df$Category[!TIMM21_Coloc_H4_Df$Category %in% names(table(TIMM21_Coloc_H4_Df$Category))[table(TIMM21_Coloc_H4_Df$Category) > 2]] <- 'Other'
# 
# ggplot2::ggplot(data = TIMM21_Coloc_H4_Df %>% dplyr::filter(PP_H4 > 0.05),aes(x=reorder(Tissue,-PP_H4),y=PP_H4)) + geom_bar(stat = 'identity') + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + facet_wrap(~Category,scales = 'free_x',nrow=1) + ylim(0,1) + ylab('Posterior - Shared Causal Variant') + xlab('Tissue') + ggtitle('Coloc - TIMM21~rs1215990')
