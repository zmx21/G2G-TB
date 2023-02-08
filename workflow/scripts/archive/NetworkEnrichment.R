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
  #Read in PASCAL files with P-values and Gene IDs
  res_list <- list()
  res_fusion_list <- list()
  lineages <- names(g2g_obj)
  for(i in 1:length(lineages)){
    cur_lineage <- lineages[i]
    cur_PASCAL <- glue::glue("{results_path}/LINEAGE_{cur_lineage}/PASCAL/{gsub(x=names(g2g_obj[[i]]),pattern='.glm.logistic.hybrid.gz',replacement='.max.genescores.txt')}")
    cur_PASCAL_fusion <- glue::glue("{results_path}/LINEAGE_{cur_lineage}/PASCAL/{gsub(x=names(g2g_obj[[i]]),pattern='.glm.logistic.hybrid.gz',replacement='.max.fusion.genescores.txt')}")
    
    for(k in 1:length(cur_PASCAL)){
      if(file.exists(cur_PASCAL[k])){
        res <- list(data.table::fread(cur_PASCAL[k]))
        res_fusion <- list(data.table::fread(cur_PASCAL_fusion[k]))
        
        mtb_gene <- strsplit(strsplit(names(g2g_obj[[i]])[k],split = ':')[[1]][1],split = '_')[[1]][2]
        names(res) <- mtb_gene
        res_list <- c(res_list,res)
        res_fusion_list <- c(res_fusion_list,res_fusion)
      }
    }
  }
  
  all_host_genes <- unique(as.character(unique(unlist(sapply(res_list,function(x) x$gene_id)))))
  #Gene-Gene Matrix
  res_matrix <- matrix(data = 0,nrow = length(res_list),ncol = length(all_host_genes))
  rownames(res_matrix) <- names(res_list)
  colnames(res_matrix) <- all_host_genes
  
  #Cluster-Cluster Matrix
  res_matrix_cluster <- matrix(data = 0,nrow = length(mtb_network_tbl$genes),ncol = length(host_gene_sets$genesets))
  rownames(res_matrix_cluster) <- names(mtb_network_tbl$genes)
  colnames(res_matrix_cluster) <- host_gene_sets$geneset.names
  
  #Binarize Associations (based on GW significant threshold)
  for(i in 1:length(res_list)){
      res_matrix[i,as.character(res_list[[i]]$gene_id)] <- ifelse(p.adjust(res_list[[i]]$pvalue,method = 'fdr') < p_thresh,1,0)
  }
  
  edges_cluster <- list()
  #Construct Cluster-Cluster Matrix
  for(p in 1:nrow(res_matrix_cluster)){
    edges_cluster[[p]] <- list()
    for(q in 1:ncol(res_matrix_cluster)){
      print(c(p,q))
    
      cur_mtb_cluster <- intersect(mtb_network_tbl$genes[[p]],rownames(res_matrix))
      cur_host_cluster <- intersect(as.character(host_gene_sets$genesets[[q]]),colnames(res_matrix))
      cur_edges <- res_matrix[cur_mtb_cluster,cur_host_cluster,drop=F]
      
      #Minimize Genes by looking at fusion (due to LD)
      cur_fusion_genes <- lapply(res_fusion_list[[p]]$gene_id,function(x) strsplit(x=x,split = '_')[[1]][2:length(strsplit(x=x,split = '_')[[1]])])
      #Remove fusion genes that don't contain cluster genes
      cur_fusion_keep <- sapply(cur_fusion_genes,function(x) length(setdiff(x,cur_host_cluster)) < 1)
      if(length(cur_fusion_keep) > 0){
        cur_fusion_genes <- cur_fusion_genes[cur_fusion_keep]
        cur_fusion_genes <- cur_fusion_genes[order(sapply(cur_fusion_genes,length),decreasing = T)]
        repeating_fusion_genes <- c()
        if(length(cur_fusion_genes) > 1){
          repeating_fusion_genes <- unlist(lapply(1:length(cur_fusion_genes),function(x) which(sapply(cur_fusion_genes[x+1:length(cur_fusion_genes)],function(y) length(intersect(cur_fusion_genes[[x]],y)) > 0)) + 1))
        }

        if(length(repeating_fusion_genes) > 0){
          cur_fusion_genes <- cur_fusion_genes[-repeating_fusion_genes] 
        }
        
        cur_edges_fusion <- do.call(cbind,lapply(cur_fusion_genes,function(x){
          tmp <- as.matrix(apply(cur_edges[,x,drop=F],1,function(x) as.numeric(any(x==1))))
          colnames(tmp) <- paste0(x,collapse = '_')
          return(tmp)
        }))
        col_to_remove <- unlist(cur_fusion_genes)
        cur_edges <- cbind(cur_edges[,!colnames(cur_edges) %in% col_to_remove,drop=F],cur_edges_fusion)
      }

      if(!any(apply(cur_edges,1,function(x) sum(x != 0)) > 1) | !any(apply(cur_edges,2,function(x) sum(x != 0)) > 1)){
        res_matrix_cluster[p,q] <- 0
      }else{
        res_matrix_cluster[p,q] <- sum(as.vector(cur_edges)) / (nrow(cur_edges) * ncol(cur_edges))
      }
      edges_cluster[[p]][[q]] <- cur_edges
    }
  }
  
  return(list(res_matrix=res_matrix,res_matrix_cluster=res_matrix_cluster,edges_cluster=edges_cluster))
}
mtb_network_tbl <- ParseMtbNetwork('~/G2G_TB/data/Mtb/MTB_Network/mtb_bicluster_attributes.txt')
host_gene_sets <- GSA.read.gmt('~/Software/PASCAL/resources/genesets/msigdb/msigBIOCARTA_KEGG_REACTOME.gmt') 
results_path <- '~/G2G_TB/results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/'
g2g_obj <- readRDS(paste0(results_path,'G2G_results.rds'))
g2g_adj_mat <- ConstructAdjMatrix(results_path,g2g_obj,mtb_network_tbl,host_gene_sets,p_thresh = 0.05)
saveRDS(g2g_adj_mat,'~/G2G_TB/results/Network/g2g_adj_mat_0p05.rds')
g2g_adj_mat <- ConstructAdjMatrix(results_path,g2g_obj,mtb_network_tbl,host_gene_sets,p_thresh = 0.1)
saveRDS(g2g_adj_mat,'~/G2G_TB/results/Network/g2g_adj_mat_0p1.rds')
g2g_adj_mat <- ConstructAdjMatrix(results_path,g2g_obj,mtb_network_tbl,host_gene_sets,p_thresh = 0.15)
saveRDS(g2g_adj_mat,'~/G2G_TB/results/Network/g2g_adj_mat_0p15.rds')
g2g_adj_mat <- ConstructAdjMatrix(results_path,g2g_obj,mtb_network_tbl,host_gene_sets,p_thresh = 0.2)
saveRDS(g2g_adj_mat,'~/G2G_TB/results/Network/g2g_adj_mat_0p2.rds')
