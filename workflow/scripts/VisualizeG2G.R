library(ggtree)
library(RColorBrewer)
library(ggplot2)
library(ape)

VisualizeG2G <- function(G2G_Obj,host_snp,AA_variant,lineage,phylotree,out_path){
  host_genotype <- snpStats::read.plink(bed = G2G_Obj$host_path[lineage],select.snps = host_snp)
  host_dosage <- as(host_genotype$genotypes, Class = 'numeric')
  host_dosage_ordered <- host_dosage[match(G2G_Obj$both_IDs_to_keep[[lineage]]$FAM_ID,rownames(host_dosage)),]
  host_dosage_ordered[host_dosage_ordered==2] <- paste0(host_genotype$map$allele.2,host_genotype$map$allele.2)
  host_dosage_ordered[host_dosage_ordered==1] <- paste0(host_genotype$map$allele.1,host_genotype$map$allele.2)
  host_dosage_ordered[host_dosage_ordered==0] <- paste0(host_genotype$map$allele.1,host_genotype$map$allele.1)
  
  pathogen_dosage <- G2G_Obj$aa_matrix_filt[[lineage]][,AA_variant,drop = F]
  pathogen_dosage[pathogen_dosage == 2] <- 1
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
      new_gg <- new_gg %<+% data.frame(metadata,Patho_Variant = factor(as.vector(pathogen_dosage_ordered[,i]))) + geom_tippoint(aes(shape = factor(Patho_Variant),color = factor(Host_SNP)),size=2,show.legend = T) + scale_color_manual(values=c('grey','green','red')) + labs(color = host_snp,shape = AA_variant[i])
      system(glue::glue("mkdir -p {out_path}/LINEAGE_{lineage}/"))
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
    # p1 = set_hilight_legend(gg, colours, names) + theme(legend.position="right")
    p1 = gg
    system(glue::glue("mkdir -p {out_path}/LINEAGE_{lineage}/"))
    ggsave(filename = glue::glue("{out_path}/LINEAGE_{lineage}/{host_snp}.pdf"),plot = p1,height = 10,width = 15)
  }
}
VisualizeG2G(G2G_Obj = readRDS('../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10//G2G_Obj.rds'),
             host_snp = 'rs12151990',
             AA_variant = 'Rv2348c_Rv2348c:2626678:p.Ile101Met',
             lineage = 'L4',
             phylotree = '../data/Mtb/RAxML_bestTree.Sinergia_final_dataset_human_bac_genome_available_rerooted.nwk',
             out_path = '../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10//PLINK/PC_3_pPC_0/Stratified_False/')

# VisualizeG2G(G2G_Obj = readRDS('../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True/G2G_Obj.rds'),
#              host_snp = 'rs75769176',
#              AA_variant = 'fixA_Rv3029c:3388671:p.Thr67Met',
#              lineage = 'L3',
#              phylotree = '../data/Mtb/RAxML_bestTree.Sinergia_final_dataset_human_bac_genome_available_rerooted.nwk',
#              out_path = '../results/Burden_False_SIFT_False_Del_False_HomoOnly_True/PLINK/PC_3_pPC_0/Stratified_True/')