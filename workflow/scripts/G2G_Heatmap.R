G2G_Res <- readRDS('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/G2G_results.rds')
G2G_Res_Simple <- lapply(G2G_Res$ALL,function(x) dplyr::select(x,'#CHROM',POS,P))

scale_factor <- 100
host_chr_length <- getChromInfoFromUCSC("hg19") %>% dplyr::filter(chrom %in% c(paste0('chr',1:22),'chrX'))
G2G_Df <- data.frame(Mtb_Pos=numeric(),Host_Pos=numeric(),P=numeric())
for(i in 1:length(G2G_Res_Simple)){
  cur_Mtb_pos <- as.numeric(strsplit(x=names(G2G_Res_Simple)[i],split = ':')[[1]][2])
  for(j in 1:nrow(G2G_Res_Simple[[i]])){
    if(nrow(G2G_Res_Simple[[i]]) == 0){
      next
    }
    cur_chrom <- G2G_Res_Simple[[i]]$`#CHROM`[j]
    if(paste0('chr',cur_chrom) == 'chr1'){
      cum_length = 0
    }else{
      cum_length <- sum(host_chr_length$size[1:(which(host_chr_length$chrom == paste0('chr',cur_chrom)) - 1)]/ scale_factor)
    }
    cur_host_pos <- G2G_Res_Simple[[i]]$POS[j]/scale_factor + cum_length
    cur_P <- G2G_Res_Simple[[i]]$P[j]
    G2G_Df <- rbind(G2G_Df,data.frame(Mtb_Pos=cur_Mtb_pos,Host_Pos=cur_host_pos,P=cur_P))
  }
}

P_Thresh <- 1.02459e-10
ggplot(data = G2G_Df) + geom_point(aes(x=Host_Pos,y=Mtb_Pos,color = -log10(P))) +
  scale_x_continuous(breaks = c(0,cumsum(host_chr_length$size / scale_factor)[1:22]),labels = host_chr_length$chrom,expand = c(0, 0),limits = c(0,max(G2G_Df$Host_Pos) + 200000)) +
  xlab('Human Genome Position') + ylab('M.tb Genome Position') + guides(shape='none')

