library(GSA)
ParseMtbNetwork <- function(tbl_path){
  #Genes of each cluster
  tbl <- data.table::fread(tbl_path)
  genes <- lapply(tbl$genes,function(x) unlist(strsplit(x=x,split = '\\|')))
  names(genes) <- tbl$bicluster
  
  #Regulators of each cluster 
  regulators <- lapply(tbl$regulator,function(x) strsplit(gsub(x=gsub(x=x,pattern = '-up-regulates-bicluster',replacement = ''),pattern = '-down-regulates-bicluster',replacement = ''),split = '\\|')[[1]])
  regulators <- lapply(regulators,function(x) {
    if(length(x) == 1){
      if(is.na(x)){
        return(NA)
      }
    }else{
      return(sapply(x,function(y) strsplit(x=y,split = '_')[[1]][1],USE.NAMES = F)) 
    }})
  names(regulators) <- tbl$bicluster
  #Info for each cluster
  tbl_info <- dplyr::select(tbl,bicluster,n_genes=nrow,score,functions)
  
  return(list(genes=genes,regulators=regulators,tbl_info=tbl_info))
}

ConstructAdjMatrix <- function(results_path,g2g_obj,mtb_network_tbl,host_gene_sets,p_thresh){
  res_list <- list()
  lineages <- names(g2g_obj)
  for(i in 1:length(lineages)){
    cur_lineage <- lineages[i]
    cur_PASCAL <- glue::glue("{results_path}/LINEAGE_{cur_lineage}/PASCAL/{gsub(x=names(g2g_obj[[i]]),pattern='.gz',replacement='.sum.genescores.txt')}")
    for(k in 1:length(cur_PASCAL)){
      if(file.exists(cur_PASCAL[k])){
        res <- list(data.table::fread(cur_PASCAL[k]))
        mtb_gene <- strsplit(strsplit(names(g2g_obj[[i]])[k],split = ':')[[1]][1],split = '_')[[1]][2]
        names(res) <- mtb_gene
        res_list <- c(res_list,res)
      }
    }
  }
  
  all_host_genes <- as.character(unique(unlist(sapply(res_list,function(x) x$gene_id))))
  res_matrix <- matrix(data = 0,nrow = length(res_list),ncol = length(all_host_genes))
  rownames(res_matrix) <- names(res_list)
  colnames(res_matrix) <- all_host_genes
  
  res_matrix_cluster <- matrix(data = 0,nrow = length(mtb_network_tbl$genes),ncol = length(host_gene_sets$genesets))
  rownames(res_matrix_cluster) <- names(mtb_network_tbl$genes)
  colnames(res_matrix_cluster) <- host_gene_sets$geneset.names
  
  for(i in 1:length(res_list)){
    res_matrix[i,as.character(res_list[[i]]$gene_id)] <- ifelse(p.adjust(res_list[[i]]$pvalue,method = 'fdr') < p_thresh,1,0)
  }
  
  for(p in 1:nrow(res_matrix_cluster)){
    print(p)
    for(q in 1:ncol(res_matrix_cluster)){
      cur_mtb_cluster <- intersect(mtb_network_tbl$genes[[p]],rownames(res_matrix))
      cur_host_cluster <- intersect(as.character(host_gene_sets$genesets[[q]]),colnames(res_matrix))
      
      cur_edges <- res_matrix[cur_mtb_cluster,cur_host_cluster,drop=F]
      if(!any(apply(cur_edges,1,function(x) sum(x != 0)) > 1) | !any(apply(cur_edges,2,function(x) sum(x != 0)) > 1)){
        res_matrix_cluster[p,q] <- 0
      }else{
        res_matrix_cluster[p,q] <- sum(as.vector(cur_edges)) / (nrow(cur_edges) * ncol(cur_edges))
      }
    }
  }
  
  return(list(res_matrix=res_matrix,res_matrix_cluster=res_matrix_cluster))
}
mtb_network_tbl <- ParseMtbNetwork('~/G2G_TB/data/Mtb/MTB_Network/mtb_bicluster_attributes.txt')
host_gene_sets <- GSA.read.gmt('~/Software/PASCAL/resources/genesets/msigdb/msigBIOCARTA_KEGG_REACTOME.gmt') 
results_path <- '~/G2G_TB/results/Burden_False_SIFT_False_Del_True_HomoOnly_True/AFGR_Tanz_ChrX/PLINK/PC_3_pPC_0/'
# g2g_obj <- readRDS(paste0(results_path,'G2G_results.rds'))
# g2g_adj_mat <- ConstructAdjMatrix(results_path,g2g_obj,mtb_network_tbl,host_gene_sets,p_thresh = 0.05)

tbscore_obj <- readRDS(paste0(results_path,'tbscore_int_results.rds'))
tbscore_adj_mat <- ConstructAdjMatrix(results_path,tbscore_obj,mtb_network_tbl,host_gene_sets,p_thresh = 0.1)

library(ggbipart)
adj_mat <- tbscore_adj_mat$res_matrix[intersect(rownames(tbscore_adj_mat$res_matrix),mtb_network_tbl$genes[[227]]),intersect(host_gene_sets$genesets[[249]],colnames(tbscore_adj_mat$res_matrix))]
g<- bip_railway(adj_mat, label=T)
