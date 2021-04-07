library(snpStats)
library(ape)
library(ggtree)
library(ggplot2)
library(scales)
library(phylolm)
library(phangorn)
library(geiger)
VisualizeG2G <- function(G2G_Obj,host_snp,AA_variant,lineage,phylotree,out_path){
  host_genotype <- snpStats::read.plink(bed = G2G_Obj$host_path[lineage],select.snps = host_snp)
  host_dosage <- as(host_genotype$genotypes, Class = 'numeric')
  host_dosage_ordered <- host_dosage[match(G2G_Obj$both_IDs_to_keep[[lineage]]$FAM_ID,rownames(host_dosage)),]
  
  pathogen_dosage <- G2G_Obj$aa_matrix_filt[[lineage]][,AA_variant,drop = F]
  pathogen_dosage_ordered <- pathogen_dosage[match(G2G_Obj$both_IDs_to_keep[[lineage]]$G_NUMBER,rownames(pathogen_dosage)),,drop = F]
  
  metadata <- data.frame(ID = rownames(pathogen_dosage_ordered),Host_SNP = host_dosage_ordered)
  rownames(metadata) <- metadata$ID
  tree <- ape::read.tree(phylotree)
  tree_filt <- ape::keep.tip(tree,rownames(pathogen_dosage_ordered))
  
  gg <- ggtree::ggtree(tree_filt,layout='circular',branch.length = 'none')
  gg <- gg %<+% metadata + geom_tippoint(aes(color = factor(Host_SNP)),size=1.5,show.legend = T) + scale_color_manual(values=c('grey','green','red')) + labs(color = host_snp)
  
  nodes <- rep(NA,length(AA_variant))
  colours <- c()
  cl_gradient <- c(brewer.pal(n = 8, name = 'Dark2'),brewer.pal(n = 11, name = 'RdYlBu'))
  names <- c()
  cl_cnt <- 1
  any_clades <- F
  for(i in 1:length(AA_variant)){
    cur_dosage <- as.vector(pathogen_dosage_ordered[,i])
    cur_minor_allele <- as.numeric(names(table(cur_dosage))[which.min(table(cur_dosage))])
    nodes[i] <- getMRCA(tree_filt, as.character(metadata$ID[which(cur_dosage==cur_minor_allele)]))
    #If AA variant is not monophyletic
    if(length(unlist(phangorn::Descendants(tree_filt,nodes[i],type = 'tips'))) != length(which(cur_dosage==cur_minor_allele))){
      new_gg <- ggtree::ggtree(tree_filt,layout='circular',branch.length = 'none')
      new_gg <- new_gg %<+% data.frame(metadata,Patho_Variant = factor(as.vector(pathogen_dosage_ordered[,i]))) + geom_tippoint(aes(shape = factor(Patho_Variant),color = factor(Host_SNP)),size=1.5,show.legend = T) + scale_color_manual(values=c('grey','green','red')) + labs(color = host_snp,shape = AA_variant[i])
      ggsave(filename = glue::glue("{out_path}/LINEAGE_{lineage}/{host_snp}_{AA_variant[i]}.pdf"),plot = new_gg,height = 10,width = 15)
      next
    }
    any_clades <- T
    if(i > 1){
      if(nodes[i] == nodes[i-1]){
        names[cl_cnt] <- paste0(names[cl_cnt],'\n',AA_variant[i])
      }else{
        cl_cnt <- cl_cnt + 1
        names <- c(names,AA_variant[i])
        gg <- gg + geom_hilight(node = nodes[i],fill = cl_gradient[cl_cnt],extend = -0.5*cl_cnt,alpha = 0.2)
        colours <- c(colours,cl_gradient[cl_cnt])
      }
    }else{
      colours <- c(colours,cl_gradient[cl_cnt])
      names <- c(names,AA_variant[i])
      gg <- gg + geom_hilight(node = nodes[i],fill = cl_gradient[cl_cnt],extend = -0.5*cl_cnt,alpha = 0.2)
    }
  }
  if(any_clades){
    p1 = set_hilight_legend(gg, colours, names) + theme(legend.position="right")
    
    ggsave(filename = glue::glue("{out_path}/LINEAGE_{lineage}/{host_snp}.pdf"),plot = p1,height = 10,width = 15)
  }
}
GetUniqResults <- function(result_df){
  do.call(rbind,lapply(1:length(result_df),function(i) {
    cur_result <- result_df[[i]]
    cbind(cur_result,
          data.frame(AA =  rep(paste0(strsplit(names(result_df)[i],split = "\\.")[[1]][1:2],collapse = '.'),nrow(cur_result)),stringsAsFactors = F))
    
  }))
}

FigtPhyloglm <- function(G2G_Obj,host_snp = NULL,AA_variant,lineage,phylotree,is_burden = T){
  if(!is.null(host_snp)){
    host_genotype <- snpStats::read.plink(bed = G2G_Obj$host_path[lineage],select.snps = host_snp)
    host_dosage <- as(host_genotype$genotypes, Class = 'numeric')
    host_dosage_ordered <- host_dosage[match(G2G_Obj$both_IDs_to_keep[[lineage]]$FAM_ID,rownames(host_dosage)),]
  }
  if(is_burden){
    pathogen_dosage <- G2G_Obj$gene_burden_non_syn[,AA_variant,drop = F]
    ids_to_keep <- rownames(pathogen_dosage)
    pathogen_dosage_ordered <- pathogen_dosage[match(G2G_Obj$both_IDs_to_keep[[lineage]]$G_NUMBER,ids_to_keep),]
  }else{
    pathogen_dosage <- G2G_Obj$aa_matrix_filt[[lineage]][,AA_variant,drop = F]
    pathogen_dosage_ordered <- pathogen_dosage[match(G2G_Obj$both_IDs_to_keep[[lineage]]$G_NUMBER,names(pathogen_dosage))]
  }

  tree <- ape::read.tree(phylotree)
  tree_filt <- ape::keep.tip(tree,ids_to_keep)
  
  if(!is.null(host_snp)){
    metadata <- data.frame(ID = names(pathogen_dosage_ordered),Host_SNP = host_dosage_ordered,Patho_Variant = pathogen_dosage_ordered)
    rownames(metadata) <- metadata$ID
    fit <- phylolm('Patho_Variant ~ Host_SNP',data = metadata,phy = tree_filt,method = 'logistic_MPLE',boot = 1000)
  }else{
    fit <- fitDiscrete(phy = tree_filt,dat = factor(pathogen_dosage_ordered),model = 'ER',transform = 'lambda')
    return(fit$opt$lambda)
  }

}
# L3_snps <- c('rs925411')

# L3_snps <- c('rs925411','rs10476842','rs55836634','rs7584562','rs868919')
# L3_hits <- lapply(L3_snps,function(x) GetUniqResults(L3_Results_AFGR) %>% dplyr::filter(ID == x))
# lapply(1:length(L3_snps),function(i) VisualizeG2G(G2G_Obj=G2G_Obj_AFGR,
#                                                             host_snp = L3_snps[i],
#                                                             AA_variant = L3_hits[[i]]$AA,
#                                                             lineage = 'L3',
#                                                             phylotree = '~/G2G_TB/data/Mtb/fasttree_TBDar_022021_1239',
#                                                             out_path <- '~/G2G_TB/results/AFGR/PLINK/PC_3_pPC_0/'))

# FigtPhyloglm(G2G_Obj_AFGR,'rs925411','fixA_Rv3029c_p.Thr67Met','L3','~/G2G_TB/data/Mtb/fasttree_TBDar_022021_1239')

# L3_Results_AFGR_Uniq <- GetUniqResults(L3_Results_AFGR)
# lapply(c(1,2,5,6,7,8),function(k) VisualizeG2G(G2G_Obj=G2G_Obj_AFGR,
#              host_snp = L3_Results_AFGR_Uniq$ID[k],
#              AA_variant = paste0(strsplit(L3_Results_AFGR_Uniq$AA[k],split = "\\.")[[1]][1:2],collapse = '.'),
#              lineage = 'L3',
#              phylotree = '~/G2G_TB/data/Mtb/fasttree_TBDar_022021_1239',
#              out_path <- '~/G2G_TB/results/AFGR/PLINK/PC_3_pPC_0/'))


# L4_Results_AFGR_Uniq <- GetUniqResults(L4_Results_AFGR)
# 
# lapply(c(1,3,4),function(k) VisualizeG2G(G2G_Obj=G2G_Obj_AFGR,
#                                                host_snp = L4_Results_AFGR_Uniq$ID[k],
#                                                AA_variant = paste0(strsplit(L4_Results_AFGR_Uniq$AA[k],split = "\\.")[[1]][1:2],collapse = '.'),
#                                                lineage = 'L4',
#                                                phylotree = '~/G2G_TB/data/Mtb/fasttree_TBDar_022021_1239',
#                                                out_path <- '~/G2G_TB/results/AFGR/PLINK/PC_3_pPC_0/'))
# 
# # L3_Results_Tanz_Uniq <- GetUniqResults(L3_Results_Tanz)
# rs925411 <- snpStats::read.plink(bed = '/home/zmxu/G2G_TB/data/WGS/WGS_Host_Data/rs925411.bed',select.snps = 'rs925411')

# burden_lambda <- pbmclapply(colnames(G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn)[x=sample(1:ncol(G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn),size = 20,replace = F)],function(x) FigtPhyloglm(G2G_Obj_AFGR_Tanz_ChrX_SIFT,AA_variant = x,lineage = 'ALL',phylotree = '~/G2G_TB/data/Mtb/fasttree_TBDar_022021_1239'),mc.cores = 10)

pPC1_p <- sapply(1:ncol(G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn),function(i) summary(glm(formula = 'y ~ PC1 + PC2 + PC3 + PC4 + PC5',data = data.frame(G2G_Obj_AFGR_Tanz_ChrX_SIFT$vir_pPCs[['ALL']],y = G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn[,i])))$coefficients['PC1','Pr(>|t|)'])
pPC2_p <- sapply(1:ncol(G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn),function(i) summary(glm(formula = 'y ~ PC1 + PC2 + PC3 + PC4 + PC5',data = data.frame(G2G_Obj_AFGR_Tanz_ChrX_SIFT$vir_pPCs[['ALL']],y = G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn[,i])))$coefficients['PC2','Pr(>|t|)'])
pPC3_p <- sapply(1:ncol(G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn),function(i) summary(glm(formula = 'y ~ PC1 + PC2 + PC3 + PC4 + PC5',data = data.frame(G2G_Obj_AFGR_Tanz_ChrX_SIFT$vir_pPCs[['ALL']],y = G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn[,i])))$coefficients['PC3','Pr(>|t|)'])
pPC4_p <- sapply(1:ncol(G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn),function(i) summary(glm(formula = 'y ~ PC1 + PC2 + PC3 + PC4 + PC5',data = data.frame(G2G_Obj_AFGR_Tanz_ChrX_SIFT$vir_pPCs[['ALL']],y = G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn[,i])))$coefficients['PC4','Pr(>|t|)'])
pPC5_p <- sapply(1:ncol(G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn),function(i) summary(glm(formula = 'y ~ PC1 + PC2 + PC3 + PC4 + PC5',data = data.frame(G2G_Obj_AFGR_Tanz_ChrX_SIFT$vir_pPCs[['ALL']],y = G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn[,i])))$coefficients['PC5','Pr(>|t|)'])

mdl <- lapply(1:ncol(G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn),function(i) glm(formula = 'y ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10',data = data.frame(G2G_Obj_AFGR_Tanz_ChrX_SIFT$vir_pPCs[['ALL']],y = G2G_Obj_AFGR_Tanz_ChrX_SIFT$gene_burden_non_syn[,i])))
