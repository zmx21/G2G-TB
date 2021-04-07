library(GSA)
library(parallel)
library(dplyr)
ParseMTBGO <- function(geneset_file,gene_ID_file){
  id_file <- data.table::fread(gene_ID_file) %>% dplyr::select(ID = V3,Gene_Name = V11)
  geneset <- GSA.read.gmt(geneset_file)
  
  geneset_GO <- sapply(geneset$geneset.descriptions,function(x) strsplit(x=x,split = 'term=')[[1]][2],USE.NAMES = F)
  geneset_df <- data.frame(GO= geneset_GO,Description = geneset$geneset.names,stringsAsFactors = F)
  
  geneset_parsed <- vector(mode = 'list',length = length(geneset$genesets))
  names(geneset_parsed) <- geneset_GO
  
  for(i in 1:length(geneset_parsed)){
    cur_gene_ID <- geneset$genesets[[i]]
    cur_gene_names <- id_file$Gene_Name[match(cur_gene_ID,id_file$ID)]
    cur_gene_names_split <- lapply(cur_gene_names,function(x) strsplit(x=x,split = '\\|')[[1]])
    cur_gene_Rv <- lapply(cur_gene_names_split,function(x){
      with_rv <- x[grepl(x=x,pattern = 'Rv') & !grepl(x=x,pattern = '\\.')]
      if(length(with_rv) == 0){
        with_rv <- NA
      }
      if(length(with_rv) > 1){
        return(with_rv[1])
      }
      return(with_rv)
    })
    geneset_parsed[[i]] <- unlist(cur_gene_Rv)
  }
  return(list(geneset = geneset_parsed,metadata = geneset_df))
}
RunPASCAL <- function(SOFTWARE_DIR,DATA_DIR,is_interaction = F){
  cur_wd <- getwd()
  pascal_dir <- glue::glue("{SOFTWARE_DIR}/PASCAL/")
  setwd(pascal_dir)
  system(glue::glue('mkdir -p {DATA_DIR}/PASCAL_OUT/'))
  system(glue::glue('ln -s {DATA_DIR}/PASCAL_OUT/ ./'))
  
  res_obj <- readRDS(glue::glue("{DATA_DIR}/results.rds"))
  if(!is_interaction){
    files_to_incl <- which(sapply(res_obj,function(x) nrow(dplyr::filter(x,ERRCODE == '.' & `#CHROM` != 'X'))) > 0)
    chr_to_incl <- lapply(res_obj[files_to_incl],function(x) unique(dplyr::filter(x,ERRCODE == '.')$`#CHROM`))
  }else{
    files_to_incl <- which(sapply(res_obj,function(x) {
      if(nrow(x) == 0){
        return(F)
      }else{
        ifelse(nrow(dplyr::filter(x,V13=='.')) == 0,F,T)
      }
    }))
    chr_to_incl <- lapply(res_obj[files_to_incl],function(x) unique(dplyr::filter(x,V13 == '.')$V1))
  }
  p_files <- sapply(names(res_obj)[files_to_incl],function(x) gsub(x=x,pattern = '.gz',replacement = '.p.gz'))
  mclapply(1:length(p_files),function(i){
    cur_file <- p_files[i]
    system(glue::glue("gunzip {DATA_DIR}/{cur_file}"))
    
    cur_chr <- chr_to_incl[[i]]
    for(chr in cur_chr){
      system(glue::glue("./Pascal --pval {DATA_DIR}/{gsub(cur_file,pattern = '.gz',replacement = '')} --customdir=resources/1kg_AFR --custom=AFR --mafcutoff=0 --genescoring=max --chr {chr} --outdir PASCAL_OUT"))
    }
    system(glue::glue("head -n 1 PASCAL_OUT/{gsub(cur_file,pattern = '.p.gz',replacement = '.max.genescores')}.chr{cur_chr[1]}.txt > PASCAL_OUT/{gsub(cur_file,pattern = '.p.gz',replacement = '.max.genescores.txt')}"))
    sapply(cur_chr,function(x) system(glue::glue("grep -v 'chromosome' PASCAL_OUT/{gsub(cur_file,pattern = '.p.gz',replacement = '.max.genescores')}.chr{x}.txt >> PASCAL_OUT/{gsub(cur_file,pattern = '.p.gz',replacement = '.max.genescores.txt')}")))
    sapply(cur_chr,function(x) system(glue::glue("rm PASCAL_OUT/{gsub(cur_file,pattern = '.p.gz',replacement = '.max.genescores')}.chr{x}.txt}")))
    
    system(glue::glue("pigz --fast {DATA_DIR}/{gsub(cur_file,pattern = '.gz',replacement = '')}"))
  },mc.cores = 5,mc.preschedule = F)
  system(glue::glue('unlink PASCAL_OUT'))
  setwd(cur_wd)
}

# MTB_GO <- ParseMTBGO('~/G2G_TB/data/Mtb/GO/MTB.geneset.txt','~/G2G_TB/data/Mtb/GO/gene_association.MTBBASE')
RunPASCAL('/home/zmxu/Software/',
          '/home/zmxu/G2G_TB/burden_interaction_results_SIFT/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_2_cov_Age-PatientSex-HIVStatus/LINEAGE_ALL/',is_interaction = T)