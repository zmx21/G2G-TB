library(BSgenome.Hsapiens.UCSC.hg19)

G2G_Results <- readRDS('../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/G2G_results.rds')$ALL
G2G_Results <- G2G_Results[sapply(G2G_Results,nrow) > 0]
Bac_Pos <- as.numeric(sapply(names(G2G_Results),function(x) strsplit(x=x,split = ':')[[1]][2]))

G2G_Df <- data.frame(Bac_Pos = numeric(),Host_Chr = numeric(),Host_Pos = numeric(),P = numeric())
chr.lengths = c(0,seqlengths(Hsapiens)[1:23])

for(i in 1:length(Bac_Pos)){
  cur_entry <- G2G_Results[[i]]
  for(j in 1:nrow(cur_entry)){
    cur_chr <- paste0('chr',cur_entry$`#CHROM`[j])
    offset <- which(names(chr.lengths) == cur_chr) - 1
    chr_offset <- sum(chr.lengths[1:offset])
    G2G_Df <- rbind(G2G_Df,data.frame(Bac_Pos = Bac_Pos[i],Host_Chr = cur_entry$`#CHROM`[j],Offset_Pos = cur_entry$POS[j] + chr_offset,Host_Pos = cur_entry$POS[j],P = cur_entry$P[j]))
  }
}

labels <- c(0,seqlengths(Hsapiens)[1:22])
names(labels) <- c(paste0('chr',seq(1,22)),'chrX')
ggplot2::ggplot(G2G_Df) + aes(x=Offset_Pos,y=Bac_Pos,color = -log10(P)) + geom_point() + scale_x_continuous(name = names(labels),breaks = cumsum(labels)) + 
  xlab('Host Position') + ylab('Mtb. Position')

G2G_Obj <- readRDS('../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds')
cor_mat <- cor(G2G_Obj$aa_matrix_full,use = 'complete.obs')^2
Variant_Pos <- sapply(rownames(cor_mat),function(x) as.numeric(strsplit(x=x,split = ':')[[1]][2]))


fixA_cor <- cor_mat["fixA_Rv3029c:3388671:p.Thr67Met",]
fixA_cor_df <- data.frame(POS = Variant_Pos,r2 = fixA_cor,Name = names(fixA_cor))

ggplot2::ggplot(fixA_cor_df) + aes(x=POS,y=r2) + geom_point() +
  ggrepel::geom_text_repel(aes(x=POS,y=r2,label = Name),data = fixA_cor_df %>% dplyr::filter(Name == 'fixA_Rv3029c:3388671:p.Thr67Met')) + xlab('Mtb. Position')


Rv2348c_cor <- cor_mat["Rv2348c_Rv2348c:2626678:p.Ile101Met",]
Rv2348c_cor_df <- data.frame(POS = Variant_Pos,r2 = Rv2348c_cor,Name = names(Rv2348c_cor))

ggplot2::ggplot(Rv2348c_cor_df) + aes(x=POS,y=r2) + geom_point() +
  ggrepel::geom_text_repel(aes(x=POS,y=r2,label = Name),data = Rv2348c_cor_df %>% dplyr::filter(Name == 'Rv2348c_Rv2348c:2626678:p.Ile101Met'))  + xlab('Mtb. Position')
