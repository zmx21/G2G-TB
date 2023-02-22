library(GSA)
library(pbmcapply)
library(igraph)
library(org.Hs.eg.db)
library(annotate)

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

ParseMtbNetworkStringDB <- function(){
  
}

ConstructAdjMatrix <- function(results_path,g2g_obj,pthresh = 0.1){
  #Read in PASCAL files with P-values and Gene IDs
  res_list <- list()
  lineages <- names(g2g_obj)
  for(i in 1:length(lineages)){
    cur_lineage <- lineages[i]
    cur_PASCAL <- glue::glue("{results_path}/LINEAGE_{cur_lineage}/PASCAL/{gsub(x=names(g2g_obj[[i]]),pattern='.glm.logistic.hybrid.gz',replacement='.max.genescores.txt')}")
    cur_PASCAL_fusion <- glue::glue("{results_path}/LINEAGE_{cur_lineage}/PASCAL/{gsub(x=names(g2g_obj[[i]]),pattern='.glm.logistic.hybrid.gz',replacement='.max.fusion.genescores.txt')}")
    
    for(k in 1:length(cur_PASCAL)){
      if(file.exists(cur_PASCAL[k])){
        res <- list(rbind(data.table::fread(cur_PASCAL[k]),data.table::fread(cur_PASCAL_fusion[k])))

        mtb_gene <- strsplit(strsplit(names(g2g_obj[[i]])[k],split = ':')[[1]][1],split = '_')[[1]][2]
        names(res) <- mtb_gene
        res_list <- c(res_list,res)
      }
    }
  }
  
  all_host_genes <- unique(c(as.character(unique(unlist(sapply(res_list,function(x) x$gene_id))))))
  #Initialize Gene-Gene Matrix
  res_matrix <- matrix(data = 0,nrow = length(res_list),ncol = length(all_host_genes))
  rownames(res_matrix) <- names(res_list)
  colnames(res_matrix) <- all_host_genes
  
  #Fill in Gene-Gene Matrix
  for(i in 1:length(res_list)){
    res_matrix[i,as.character(res_list[[i]]$gene_id)] <- res_list[[i]]$pvalue
  }
  res_matrix_adj <- matrix(p.adjust(res_matrix,method = 'fdr'),ncol = ncol(res_matrix))
  res_matrix_binary <- matrix(data = NA, nrow = nrow(res_matrix_adj),ncol = ncol(res_matrix_adj))
  res_matrix_binary[res_matrix_adj < pthresh] <- 1
  res_matrix_binary[res_matrix_adj >= pthresh ] <- 0
  rownames(res_matrix_binary) <- rownames(res_matrix)
  colnames(res_matrix_binary) <- colnames(res_matrix)
  
  return(list(res_matrix=res_matrix,res_matrix_binary = res_matrix_binary))
}
#Merge genes in clusters to metagenes
ConstructMetageneClusters <- function(host_gene_sets,host_genes){
  host_gene_sets_meta <- list()
  fusion_genes <- lapply(host_genes[sapply(host_genes,function(x) grepl(x=x,pattern = '_'))],function(x) strsplit(x=x,split = '_')[[1]][2:length(strsplit(x=x,split = '_')[[1]])])
  for(i in 1:length(host_gene_sets$genesets)){
    cur_host_cluster <- unlist(host_gene_sets$genesets[i])
    #Remove fusion genes that don't contain cluster genes
    cur_fusion_keep <- sapply(fusion_genes,function(x) length(intersect(x,cur_host_cluster)) > 1)
    if(sum(cur_fusion_keep) > 0){
      cur_fusion_genes <- fusion_genes[cur_fusion_keep]
      cur_fusion_genes <- cur_fusion_genes[order(sapply(cur_fusion_genes,function(x) length(intersect(x,cur_host_cluster))),decreasing = T)]
      repeating_fusion_genes <- c()
      #Find fusion genes with less overlap with current cluster
      if(length(cur_fusion_genes) > 1){
        repeating_fusion_genes <- lapply(1:length(cur_fusion_genes),function(x) which(sapply(cur_fusion_genes[x+1:length(cur_fusion_genes)],function(y) length(intersect(cur_fusion_genes[[x]],y)) > 0)) + x)
      }
      #Remove fusion genes with less overlap
      if(length(unique(unlist(repeating_fusion_genes))) > 0){
        cur_fusion_genes <- cur_fusion_genes[-unique(unlist(repeating_fusion_genes))] 
      }
      #Combine fusion genes (with maximal overlap) and single genes
      genes_to_remove <- intersect(cur_host_cluster,unlist(cur_fusion_genes))
      host_gene_sets_meta[[i]] <- c(sapply(cur_fusion_genes,function(x) paste0(x=x,collapse = '_')),setdiff(cur_host_cluster,genes_to_remove))
      
    }else{
      host_gene_sets_meta[[i]] <- cur_host_cluster
    }
  }
  names(host_gene_sets_meta) <- host_gene_sets$geneset.names
  return(host_gene_sets_meta)
  
}
ConstructClusterClusterMatrix <- function(g2g_adj_mat,host_gene_sets_meta,mtb_network_gene_sets,store_edges = T){
  res_matrix_cluster <- matrix(NA,nrow = length(mtb_network_gene_sets),ncol = length(host_gene_sets_meta))
  rownames(res_matrix_cluster) <- names(mtb_network_gene_sets)
  colnames(res_matrix_cluster) <- names(host_gene_sets_meta)
  edges_cluster <- list()
  #Construct Cluster-Cluster Matrix
  for(p in 1:nrow(res_matrix_cluster)){
    #print(p)
    cur_mtb_cluster <- intersect(mtb_network_gene_sets[[p]],rownames(g2g_adj_mat))
    #Loop through all host clusters
    cur_res <- lapply(1:ncol(res_matrix_cluster),function(q){
      cur_host_cluster <- intersect(host_gene_sets_meta[[q]],colnames(g2g_adj_mat))
      cur_edges <- g2g_adj_mat[cur_mtb_cluster,cur_host_cluster,drop=F]
      
      if(!any(apply(cur_edges,1,function(x) sum(x != 0)) > 1) | !any(apply(cur_edges,2,function(x) sum(x != 0)) > 1)){
        cur_density <- 0
      }else{
        cur_density <- sum(as.vector(cur_edges)) / (nrow(cur_edges) * ncol(cur_edges))
      }
        return(list(cur_edges=cur_edges,cur_density = cur_density))
    })
    if(store_edges){
      edges_cluster[[p]] <- lapply(cur_res,function(x) x$cur_edges)
    }
    res_matrix_cluster[p,] <- sapply(cur_res,function(x) x$cur_density)
  }
  if(store_edges){
    return(list(res_matrix_cluster=res_matrix_cluster,edges_cluster=edges_cluster))
  }else{
    return(list(res_matrix_cluster=res_matrix_cluster))
  }
}

PermuteGeneset <- function(seed,input_geneset){
  set.seed(seed)
  all_genes <- unlist(input_geneset)
  cluster_sizes <- sapply(input_geneset,length)
  
  #permute overall gene order
  permuted_genes <- all_genes[sample(1:length(all_genes),size = length(all_genes),replace = F)]
  
  #construct new permuted geneset
  output_geneset <- list()
  cnt = 1
  for(i in 1:length(input_geneset)){
    cur_length <- cluster_sizes[i]
    cur_perm_genes <- permuted_genes[cnt:(cnt + cur_length - 1)]
    names(cur_perm_genes) <- NULL
    output_geneset[[i]] <- cur_perm_genes
    cnt <- cnt + cur_length
  }
  
  names(output_geneset) <- paste0('PERM',seed,'.',names(input_geneset))
  return(output_geneset)
}

VisualizeClusters <- function(edge_matrix){
  #Convert ENTREZ to Symbol
  colnames(edge_matrix) <- unlist(lookUp(colnames(edge_matrix), 'org.Hs.eg', 'SYMBOL'))
  
  edges <- which(edge_matrix == 1,arr.ind = T)
  Mtb_Nodes <- rownames(edge_matrix)[edges[,'row',drop=T]]
  Host_Nodes <- colnames(edge_matrix)[edges[,'col',drop=T]]
  edges_df <- data.frame(Mtb_Nodes=Mtb_Nodes,Host_Nodes=Host_Nodes)
  
  g <- graph.data.frame(edges_df, directed=FALSE,vertices = data.frame(c(rownames(edge_matrix),colnames(edge_matrix))))
  V(g)$type <- c(rep(T,nrow(edge_matrix)),rep(F,ncol(edge_matrix)))  
  
  V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
  V(g)$shape <- ifelse(V(g)$type, "circle", "square")
  
  
  plot(g, layout=layout.bipartite, vertex.size=20, vertex.label.cex=0.8)
  
}



# results_path <- '~/G2G_TB/results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/'
# g2g_obj <- readRDS(paste0(results_path,'G2G_results.rds'))
# g2g_adj_mat <- ConstructAdjMatrix(results_path,g2g_obj,pthresh = 0.15)
# saveRDS(g2g_adj_mat,'~/G2G_TB/results/Network/G2G_Adj_Mat.rds')
# 
# host_gene_sets <- GSA.read.gmt('~/Software/PASCAL/resources/genesets/msigdb/msigBIOCARTA_KEGG_REACTOME.gmt')
# host_gene_sets_meta <- ConstructMetageneClusters(host_gene_sets,colnames(g2g_adj_mat$res_matrix_binary))
# saveRDS(host_gene_sets_meta,'~/G2G_TB/results/Network/msigBIOCARTA_KEGG_REACTOME_meta.rds')
# 
# mtb_network_tbl <- ParseMtbNetwork('~/G2G_TB/data/Mtb/MTB_Network/mtb_bicluster_attributes.txt')
# mtb_network_gene_sets <- mtb_network_tbl$genes
# g2g_cluster_mat <- ConstructClusterClusterMatrix(g2g_adj_mat = g2g_adj_mat$res_matrix_binary,
#                                                  host_gene_sets_meta = host_gene_sets_meta,
#                                                  mtb_network_gene_sets = mtb_network_gene_sets)
# saveRDS(g2g_cluster_mat,'~/G2G_TB/results/Network/G2G_Cluster_Mat.rds')
# 
# g2g_adj_mat <- readRDS('~/G2G_TB/results/Network/G2G_Adj_Mat.rds')
# host_gene_sets_meta <- readRDS('~/G2G_TB/results/Network/msigBIOCARTA_KEGG_REACTOME_meta.rds')
# mtb_network_tbl <- ParseMtbNetwork('~/G2G_TB/data/Mtb/MTB_Network/mtb_bicluster_attributes.txt')
# mtb_network_gene_sets <- mtb_network_tbl$genes
# 
# N_Perm <- 1000
# N_Cores <- 10
# permuted_host_geneset <- lapply(1:N_Perm,function(i) PermuteGeneset(i,host_gene_sets_meta))
# saveRDS(permuted_host_geneset,'~/G2G_TB/results/Network/perm_host.rds')
# permuted_mtb_geneset <- lapply(1:N_Perm,function(i) PermuteGeneset(i,mtb_network_gene_sets))
# saveRDS(permuted_mtb_geneset,'~/G2G_TB/results/Network/perm_mtb.rds')
# 
# g2g_cluster_mat_perm <- pbmclapply(1:N_Perm,function(i) ConstructClusterClusterMatrix(g2g_adj_mat = g2g_adj_mat$res_matrix_binary,
#                                                                                     host_gene_sets_meta = permuted_host_geneset[[i]],
#                                                                                     mtb_network_gene_sets = permuted_mtb_geneset[[i]],
#                                                                                     store_edges = F),mc.cores = N_Cores)
# saveRDS(g2g_cluster_mat_perm,'~/G2G_TB/results/Network/perm_results.rds')
# 
# g2g_cluster_mat <- readRDS('~/G2G_TB/results/Network/G2G_Cluster_Mat.rds')
# VisualizeClusters(g2g_cluster_mat$edges_cluster[[582]][[7]])
# 
# 
# g2g_cluster_mat_perm <- readRDS('~/G2G_TB/results/Network/perm_results.rds')
# perm_density <- sapply(g2g_cluster_mat_perm,function(x) max(x$res_matrix_cluster))

# Density estimation
den <- density(perm_density)

# Plot
hist(perm_density,freq = FALSE, col="gray", breaks = 20,main = 'Permutation Density Values',xlab = 'Module-Module Density')
lines(den)

value <- max(g2g_cluster_mat$res_matrix_cluster)
polygon(c(den$x[den$x >= value ], value),
        c(den$y[den$x >= value ], 0),
        col = "red",
        border = 1)
text(x = 0.15,y = 5,paste0('p = ',sum(perm_density > value) / length(perm_density)))