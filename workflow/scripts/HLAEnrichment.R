library(seqinr)
library(Biostrings)
library(UniprotR)
library(dplyr)
library(ggpubr)
library(Matching)

GetResults <- function(OUT_DIR,suffix = 'glm.logistic.hybrid',p_thresh=5e-8,n_cores=5,is_interaction = F,is_ordinal = F,tool = 'PLINK'){
  all_files <- dir(OUT_DIR)
  result_files <- all_files[sapply(all_files,function(x) grepl(pattern = suffix,x=x))]
  if(!is_interaction){
    if((tool == 'PLINK' | tool == 'PLINK-FIRTH') & !grepl(x=suffix,pattern = 'gz')){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0("awk '{ if (NR == 1 || $13 <= ",p_thresh,") {print} }' ",OUT_DIR,x)),mc.cores = n_cores)
    }else if ((tool == 'PLINK' | tool == 'PLINK-FIRTH') & grepl(x=suffix,pattern = 'gz')){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0('zcat ',OUT_DIR,x," | awk '{ if (NR == 1 || $13 <= ",p_thresh,") {print} }' ")),mc.cores = n_cores)
    }else if(tool == 'HLA-PLINK' ){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0('zcat ',OUT_DIR,x)),mc.cores = n_cores)
    }else if(tool == 'HLA-PLINK-PERM' ){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0('zcat ',OUT_DIR,x),select = c('ID','P')),mc.cores = n_cores)
    }
  }else{
    if(is_ordinal){
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0("zcat ",OUT_DIR,x," | awk -F ',' '{ if ($6 <= ",p_thresh,") {print} }' ")),mc.cores = n_cores)
    }else{
      results <- pbmclapply(result_files,function(x) data.table::fread(cmd = paste0("zcat ",OUT_DIR,x," | grep 'ADDx' | awk '{ if ($12 <= ",p_thresh,") {print} }' ")),mc.cores = n_cores)
    }
  }
  names(results) <- result_files
  return(results)
}

ParseHLAResult <- function(HLA_Results, type){
  if(type == 'AA'){
    HLA_Results$HLA_Pos <- sapply(HLA_Results$ID,function(x) paste0(strsplit(x=x,split = '_')[[1]][1:5],collapse = '_'))
    single_mutation_only <- sapply(HLA_Results$ID,function(y) nchar(strsplit(x = y,split = '_')[[1]][6]) == 1)
    HLA_Results <- HLA_Results[single_mutation_only,]
    results_df <- dplyr::group_by(HLA_Results,HLA_Pos)  %>% 
      filter(P == min(P)) %>% 
      dplyr::arrange(P) %>%
      filter(1:n() == 1) %>% dplyr::select(HLA_Pos,OR,P)
  }else if (type == 'Allele'){
    digits <- as.numeric(sapply(HLA_Results$ID,function(x) length(strsplit(x=x,split = ':')[[1]])))
    HLA_Results$digits <- digits
    results_df <- HLA_Results %>% dplyr::filter(digits == 2) #dplyr::filter(digits <= 2)
    results_df <- results_df %>% dplyr::select(Allele=ID,OR=OR,P=P)
  }
  return(results_df)
}
#Group perfectly linked variants together
ConstructVariantSets <- function(pheno_mat){
  #Correlation matrix of the phenotypes
  cor_mat <- cor(pheno_mat,use = 'complete.obs')
  
  locus_to_test <- colnames(pheno_mat)
  variant_sets <- list()
  counter = 1
  while(length(locus_to_test) > 0){
    cur_locus <- locus_to_test[1]
    cor_with_locus <- cor_mat[cur_locus,]^2
    linked_variants <- names(cor_with_locus)[cor_with_locus == 1]
    variant_sets[[counter]] <- linked_variants
    locus_to_test <- setdiff(locus_to_test[-1],linked_variants)
    counter <- counter + 1
  }
  return(variant_sets)
}

ParseEpitope <- function(FASTA_List,Epitope_File){
  results <- data.frame(Gene = character(),Start = numeric(),End = numeric(),Peptide = character())
  for(i in 1:nrow(Epitope_File)){
    print(i)
    if(Epitope_File$Gene_ID[i] %in% names(Mtb_FASTA_List)){
      cur_alignment <- as.matrix(pairwiseAlignment(Epitope_File$Peptide[i],Mtb_FASTA_List[[Epitope_File$Gene[i]]],type = 'overlap'))
      if(!all(cur_alignment == '-')){
        start <- min(which(cur_alignment != '-'))
        end <- max(which(cur_alignment != '-'))
        if((end - start + 1) == length(strsplit(Epitope_File$Peptide[i],split = '')[[1]])){
          results <- rbind(results,data.frame(Gene = Epitope_File$Gene[i],Start = start,End = end,Peptide = Epitope_File$Peptide[i]))
        }
      }
    }
  }
  return(results)
}

Epitope_File <- data.table::fread('../../data/Mtb/IEDB_Epitopes_25_10_2022.csv') %>%
  dplyr::filter(`Organism Name` == 'Mycobacterium tuberculosis H37Rv' | `Organism Name` == 'Mycobacterium tuberculosis') %>%
  dplyr::filter(`Parent Protein Accession` != '') %>% dplyr::select(ID = `Epitope ID`,Peptide = `Description`,Start = `Starting Position`,
                                                                    End = `Ending Position`,Antigen_Name = `Antigen Name`,
                                                                    Protein = `Parent Protein`,Accession = `Parent Protein Accession`) %>%
  dplyr::filter(!grepl('PE',Antigen_Name) & !grepl('PPE',Antigen_Name) & !grepl('PPE',Protein) & !grepl('PE',Protein) & !grepl('transposase',Antigen_Name) & !grepl('transposase',Protein) & !grepl("[^A-Za-z]",Peptide))
write(Epitope_File$Accession,file = '../../data/Mtb/IEDB_Epitopes_25_10_2022_Accession.txt')
Accession_Mapping <- data.table::fread('../../data/Mtb/IEDB_Epitopes_25_10_2022_Accession_Mapping.tsv') %>% dplyr::select(Accession = From,Gene_ID = To)
Epitope_File <- Epitope_File %>% dplyr::left_join(Accession_Mapping)

Mtb_FASTA <- seqinr::read.fasta('../../data/Mtb/H37Rv_AA.txt')
Mtb_FASTA_Names <- sapply(Mtb_FASTA,function(x) getAnnot(x))
Mtb_FASTA_Names <- sapply(Mtb_FASTA_Names,function(x) strsplit(x=strsplit(x=x,split = 'gene=')[[1]][2],split = ']')[[1]][1])
Mtb_FASTA_List <- lapply(Mtb_FASTA,function(x) paste0(toupper(x),collapse = ''))
names(Mtb_FASTA_List) <- Mtb_FASTA_Names
Epitope_File_Parsed <- ParseEpitope(Mtb_FASTA_List,Epitope_File)

G2G_HLA_Results <- readRDS('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/HLA-PLINK/PC_3_pPC_0/Stratified_False/G2G_results.rds')

G2G_HLA_AA <- lapply(G2G_HLA_Results$HLA_AA$ALL,function(x) ParseHLAResult(HLA_Results = x,type = 'AA'))
names(G2G_HLA_AA) <- names(G2G_HLA_Results$HLA_AA$ALL)
G2G_HLA_Allele <- lapply(G2G_HLA_Results$HLA_Allele$ALL,function(x) ParseHLAResult(HLA_Results = x,type = 'Allele'))
names(G2G_HLA_Allele) <- names(G2G_HLA_Results$HLA_Allele$ALL)

G2G_HLA_Results_Genes <- sapply(names(G2G_HLA_AA),function(x) strsplit(x=x,split = '_')[[1]][1],USE.NAMES = F)
G2G_HLA_Results_Mutations <- sapply(names(G2G_HLA_AA),function(x) strsplit(strsplit(x=x,split = ':')[[1]][3],split = '\\.')[[1]][2],USE.NAMES = F) 
G2G_HLA_Results_Pos <- as.numeric(sapply(names(G2G_HLA_AA),function(x) gsub("[^0-9.-]", "", strsplit(strsplit(x=x,split = ':')[[1]][3],split = '\\.')[[1]][2]),USE.NAMES = F))
G2G_HLA_Results_Is_Epitope <- sapply(1:length(G2G_HLA_Results_Genes),function(i){
  is_epitope <- F
  if(G2G_HLA_Results_Genes[i] %in% Epitope_File_Parsed$Gene){
    gene_index <- which(Epitope_File_Parsed$Gene == G2G_HLA_Results_Genes[i] )
    for(k in 1:length(gene_index)){
      if(G2G_HLA_Results_Pos[i] >= Epitope_File_Parsed$Start[gene_index[k]] & G2G_HLA_Results_Pos[i] <= Epitope_File_Parsed$End[gene_index[k]]){
        is_epitope <- T
      }
    }
  }
  return(is_epitope)
  
})
G2G_HLA_Results_Epitope_Index <- lapply(1:length(G2G_HLA_Results_Genes),function(i){
  epitope_index <- c()
  if(G2G_HLA_Results_Genes[i] %in% Epitope_File_Parsed$Gene){
    gene_index <- which(Epitope_File_Parsed$Gene == G2G_HLA_Results_Genes[i] )
    for(k in 1:length(gene_index)){
      if(G2G_HLA_Results_Pos[i] >= Epitope_File_Parsed$Start[gene_index[k]] & G2G_HLA_Results_Pos[i] <= Epitope_File_Parsed$End[gene_index[k]]){
        epitope_index <- c(epitope_index,gene_index[k])
      }
    }
  }
  return(epitope_index)
  
})



G2G_HLA_Results_Summary <- data.frame(Gene = G2G_HLA_Results_Genes,
                                      Position = G2G_HLA_Results_Pos,
                                      Mutation = G2G_HLA_Results_Mutations,
                                      Is_Epitope = G2G_HLA_Results_Is_Epitope,
                                      Top_P_AA = sapply(G2G_HLA_AA,function(x) min(x$P)),
                                      SD_P_AA = sapply(G2G_HLA_AA,function(x) sd(x$P) / mean(x$P)),
                                      Top_P_Allele = sapply(G2G_HLA_Allele,function(x) min(x$P)),
                                      Top_OR_AA = sapply(G2G_HLA_AA,function(x) x$OR[which.min(x$P)]),
                                      Top_OR_Allele = sapply(G2G_HLA_Allele,function(x) x$OR[which.min(x$P)]),
                                      Top_AA = sapply(G2G_HLA_AA,function(x) x$HLA_Pos[which.min(x$P)]),
                                      Top_Allele = sapply(G2G_HLA_Allele,function(x) x$Allele[which.min(x$P)]),
                                      Epitope = sapply(G2G_HLA_Results_Epitope_Index,function(x) paste0(Epitope_File_Parsed$Peptide[x],collapse = ','),USE.NAMES = F))
G2G_Obj <- readRDS('../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds')
Variant_Set <- ConstructVariantSets(G2G_Obj$aa_matrix_full)

Grouped_HLA_Results <- list()
for(i in 1:length(Variant_Set)){
  cur_variant_set <- Variant_Set[[i]]
  cur_rownames <- sapply(cur_variant_set,function(x) paste0(x,'.',x,'.glm.logistic.hybrid.gz'))
  cur_results <- G2G_HLA_Results_Summary[cur_rownames,]
  Grouped_HLA_Results[[i]] <- cur_results %>% dplyr::mutate(Gene_Set = paste0(cur_variant_set,collapse = ','))
}

Grouped_HLA_Results_Non_Epitope <- Grouped_HLA_Results[sapply(Grouped_HLA_Results,function(x) !any(x$Is_Epitope))]
Grouped_HLA_Results_Epitope <- Grouped_HLA_Results[sapply(Grouped_HLA_Results,function(x) any(x$Is_Epitope))]

Grouped_HLA_Results_Non_Epitope_df <- do.call(rbind,lapply(Grouped_HLA_Results_Non_Epitope,function(x) x %>% dplyr::select(-Position,-Mutation,-Is_Epitope,-Epitope) %>% dplyr::mutate(Epitope_Variant_Set = 'Not Annotated as \n T cell Epitope')))
Grouped_HLA_Results_Epitope_df <- do.call(rbind,lapply(Grouped_HLA_Results_Epitope,function(x) x %>% dplyr::filter(Gene %in% G2G_HLA_Results_Summary$Gene[G2G_HLA_Results_Summary$Is_Epitope] & Position %in% G2G_HLA_Results_Summary$Position[G2G_HLA_Results_Summary$Is_Epitope]) %>% 
                                                         dplyr::select(-Position,-Mutation,-Is_Epitope,-Epitope) %>% dplyr::mutate(Epitope_Variant_Set = 'T cell Epitope')))

Grouped_HLA_Results_Epitope_df$Is_Antigen <- Grouped_HLA_Results_Epitope_df$Gene %in% Epitope_File$Gene_ID
Grouped_HLA_Results_Non_Epitope_df$Is_Antigen <- Grouped_HLA_Results_Non_Epitope_df$Gene %in% Epitope_File$Gene_ID


ggplot2::ggplot(rbind(Grouped_HLA_Results_Non_Epitope_df,Grouped_HLA_Results_Epitope_df),aes(x=Epitope_Variant_Set,y=-log10(Top_P_AA)))  + geom_dotplot(binaxis='y', stackdir='centerwhole',dotsize = 0.1,color = 'red',fill = 'red') +
  ggrepel::geom_label_repel(aes(x=Epitope_Variant_Set,y=-log10(Top_P_AA),label = Gene),size = 3,data = Grouped_HLA_Results_Epitope_df,max.overlaps = 20) + xlab('Annotation') + ylab('-log10(P)') + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
                 geom="pointrange", color="black") + 
  stat_compare_means()

# ggplot2::ggplot(rbind(Grouped_HLA_Results_Non_Epitope_df,Grouped_HLA_Results_Epitope_df),aes(x=Epitope_Variant_Set,y=-log10(Top_P_Allele)))  + geom_dotplot(binaxis='y', stackdir='centerwhole',dotsize = 0.1,color = 'red') +
#   ggrepel::geom_label_repel(aes(x=Epitope_Variant_Set,y=-log10(Top_P_Allele),label = Gene),size = 3,data = Grouped_HLA_Results_Epitope_df,max.overlaps = 20) + xlab('Annotation') + ylab('-log10(P)') + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
#                                                                                                                                                                                                                  geom="pointrange", color="black") + 
#   stat_compare_means()


# N_Perm <- 100
# Epitope_Variants <- rownames(G2G_HLA_Results_Summary)[G2G_HLA_Results_Summary$Is_Epitope == T]
# HLA_Out_Path <- '../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/HLA-PLINK-PERM/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/'
# min_P_AA_Perm <- c()
# min_P_Allele_Perm <- c()
# 
# for (k in 1:N_Perm){
#   results_allele <- GetResults(paste0(HLA_Out_Path,glue::glue("HLA_Allele/Perm_{k}/")),suffix = 'glm.logistic.hybrid.gz',p_thresh=1,n_cores=20,is_interaction = F,is_ordinal = F,tool = 'HLA-PLINK-PERM')
#   results_AA <- GetResults(paste0(HLA_Out_Path,glue::glue("HLA_AA/Perm_{k}/")),suffix = 'glm.logistic.hybrid.gz',p_thresh=1,n_cores=20,is_interaction = F,is_ordinal = F,tool = 'HLA-PLINK-PERM')
#   
#   min_P_AA_Perm <- c(min_P_AA_Perm,min(sapply(results_AA[Epitope_Variants],function(x) {
#     df <- x[sapply(x$ID,function(y) nchar(strsplit(x = y,split = '_')[[1]][6]) == 1),]
#     return(min(df$P))
#     })))
#   min_P_Allele_Perm <- c(min_P_Allele_Perm,min(sapply(results_allele[Epitope_Variants],function(x) {
#     return(min(x$P))
#   })))
# }



# G2G_HLA_Results_Summary_Epitope <- G2G_HLA_Results_Summary %>% dplyr::filter(Is_Epitope)
# rownames(G2G_HLA_Results_Summary_Epitope) <- NULL
# 
# ggplot(G2G_HLA_Results_Summary,aes(x=Is_Epitope,y=-log10(Top_P_AA))) + geom_boxplot() 
# 
# HLA_AA_P <- data.frame(AA_P = unlist(sapply(G2G_HLA_AA,function(x) x$P)),
#                        Is_Epitope = unlist(sapply(1:nrow(G2G_HLA_Results_Summary),function (i) rep(x=G2G_HLA_Results_Summary$Is_Epitope[i],nrow(G2G_HLA_AA[[i]])))),
#                        Is_Min = unlist(sapply(1:nrow(G2G_HLA_Results_Summary),function (i) G2G_HLA_AA[[i]]$P == min(G2G_HLA_AA[[i]]$P))))
# 
# 
# arrow_pos <- data.frame(name = G2G_HLA_Results_Summary_Epitope$Gene,
#                         x = -log10(G2G_HLA_Results_Summary_Epitope$Top_P_AA),
#                         y = 2000)
# ggplot(HLA_AA_P)+ geom_histogram(aes(x = -log10(AA_P)),binwidth = 0.02) +
#   ggrepel::geom_text_repel(data=arrow_pos, aes(label=name, x=x, y=y), size=3,direction = 'y',min.segment.length	 = 10) + ggtitle('G2G - HLA AA and M.Tb Variants') + 
#   ylab('count')
# 
# HLA_Allele_P <- data.frame(Allele_P = unlist(sapply(G2G_HLA_Allele,function(x) x$P)))
# arrow_pos <- data.frame(name = G2G_HLA_Results_Summary_Epitope$Gene,
#                         x = -log10(G2G_HLA_Results_Summary_Epitope$Top_P_Allele),
#                         y = 2000)
# ggplot(HLA_Allele_P)+ geom_histogram(aes(x = -log10(Allele_P)),binwidth = 0.02) +
#   ggrepel::geom_text_repel(data=arrow_pos, aes(label=name, x=x, y=y), size=3,direction = 'y',min.segment.length	 = 10) + ggtitle('G2G - HLA Allele and M.Tb Variants')+ 
#   ylab('count')
# 
# 
# 
# P_Vals <- lapply(G2G_HLA_AA,function(x) x$P)

