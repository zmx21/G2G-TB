library(STRINGdb)
library(GeneTonic)
library(parallel)
library(dplyr)

#Download M.tb STRINGdb, default score threhsold (400)
string_db <- STRINGdb$new(species=83332, input_directory="../../data/",version="11.5",score_threshold=400) 
#Extract StringDB ID to gene ID mapping
string_id_df <- data.table::fread('../../data/Mtb/MTB_Network/StringDB/83332.protein.info.v11.5.txt') %>%
  dplyr::select(STRING_ID=`#string_protein_id`,Gene_ID = preferred_name)
string_id_df$Rv_ID <- sapply(string_id_df$STRING_ID,function(x) strsplit(x=x,split = '\\.')[[1]][2])

#Remove PE,PPE,PGRS, and phages
locus_to_excl <- data.table::fread('../../data/Mtb/Locus_to_exclude_Mtb.txt',fill = T)
locus_to_excl_RvIDs <- sapply(locus_to_excl$`locus tag`,function(x) ifelse(grepl(x=x,pattern = '-'),strsplit(x=x,split = '-')[[1]][2],x),USE.NAMES = F)
locus_to_excl_RvIDs <- sapply(locus_to_excl_RvIDs,function(x) ifelse(grepl(x=x,pattern = 'PE'),strsplit(x=x,split = 'PE')[[1]][1],x),USE.NAMES = F)

string_id_df <- string_id_df %>% dplyr::filter(!Rv_ID %in% locus_to_excl_RvIDs)

#Returns igraph object containing network of specified genes
GetSubNetwork <- function(string_db,STRING_IDs){
  #Create graph with specified genes
  string_subgraph <- string_db$get_subnetwork(STRING_IDs)
  string_subgraph <- igraph::simplify(string_subgraph)
}

#Loops through specified inflation parameters, returns clustering with most clusters within size limit
ChooseInflation <- function(cur_subnetwork,mcl_inflation_param = c(1.5,2,2.5,5,10,15,20,30),size_limits = c(5,30)){
  #Return number of clusters within specified size limit
  GetNumClustersWithinSize <- function(string_obj){
    n_in_limit <- sum(table(string_obj$membership) >= size_limits[1] & table(string_obj$membership) <= size_limits[2])
    return(n_in_limit)
  }
  
  #Run MCL across all specified inflation parameters
  mcl_results <- mclapply(mcl_inflation_param,function(x) GeneTonic::cluster_markov(cur_subnetwork, allow_singletons = T, mcl_inflation = x),mc.cores = length(mcl_inflation_param))
  
  #Choose inflation with max number of clusters within size limits
  optimal_inflation_ind <- which.max(sapply(mcl_results,function(x) GetNumClustersWithinSize(x)))
  print(paste0('tested inflation:',paste0(mcl_inflation_param,collapse = ','),' best = ',mcl_inflation_param[optimal_inflation_ind]))
  optimal_mcl_clusters <- mcl_results[[optimal_inflation_ind]]
  
  #Convert cluster object to list
  optimal_mcl_clusters_list <- lapply(unique(optimal_mcl_clusters$membership),function(x) optimal_mcl_clusters$names[optimal_mcl_clusters$membership == x])
  
  return(optimal_mcl_clusters_list)
  
}

## Runs MCL algorithm recursively. Stops when
#1. no improvement of clusters that fall within size limit or all clusters are within size limit OR
#2. when number of iterations exceeds limit 
#@Input
#size_limits: minimum and maximum cluster size to be considered biologically "plausible". Only clusters within this limit are retained
#string_db: downloaded stringDB database
#genes_to_incl: genes to include in network
#iter_limit: stop recursion when this is reached
RunMCL <- function(string_db,genes_to_incl,iter_limit = 2,cur_iter = 1,size_limits = c(5,30)){
  #Run MCL recursively
  if(cur_iter <= iter_limit){
    #Call MCL, loop through specified inflation values. Choose value that produce most number of clusters within limit
    cur_subnetwork <- GetSubNetwork(string_db,genes_to_incl)
    cur_clusters_list <- ChooseInflation(cur_subnetwork,size_limits = size_limits)
    
    #If MCL clustering doesn't change anything, drop cluster
    if(length(cur_clusters_list) == 1){
      print('MCL failed to Re-cluster')
      return(list())
    }
    cur_clusters_sizes <- sapply(cur_clusters_list,length)
    
    #Find clusters that exceed limit
    overlimit_clusters <- cur_clusters_list[cur_clusters_sizes > size_limits[2]]
    print(paste0(length(overlimit_clusters),' overlimit clusters'))
    
    #Find clusters that are within limit
    within_limit_clusters <- cur_clusters_list[cur_clusters_sizes <= size_limits[2] & cur_clusters_sizes >= size_limits[1]]
    print(paste0(length(within_limit_clusters),' within-limit clusters'))
    
    if(length(overlimit_clusters) > 0){
      #Run MCL for all clusters that are overlimit, appending it to the list of withinlimit clusters.
      return(c(within_limit_clusters,unlist(lapply(overlimit_clusters,function(x) RunMCL(string_db = string_db,genes_to_incl = x,
                                                                                           iter_limit = iter_limit,cur_iter = cur_iter + 1,size_limits=size_limits)),recursive = F)))
    }#If no overlimit clusters, return all clusters that are within limit. 
    else{
      return(within_limit_clusters)
    }

  }
  #Stop code if iteration exceeds limit. Does not return the oversized cluster
  else{
    print(paste0('Iteration Limit Reached'))
    return(list())
  }
}
mcl_clusters <- RunMCL(string_db = string_db,genes_to_incl = string_id_df$STRING_ID,iter_limit = 4)
#Change StringDB ID to Gene ID
mcl_clusters_recode <- lapply(mcl_clusters,function(x) string_id_df$Gene_ID[match(x,string_id_df$STRING_ID)])
#Save geneset files
saveRDS(mcl_clusters,file = '../../data/Mtb/MTB_Network/StringDB/mcl_clusters.rds')
saveRDS(mcl_clusters_recode,file = '../../data/Mtb/MTB_Network/StringDB/mcl_clusters_recode.rds')
