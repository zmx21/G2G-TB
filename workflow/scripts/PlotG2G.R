library(dplyr)
library(GenomicFeatures)
library(ggplot2)
library(latex2exp)
library(ggrepel)
library(qqman)

#Pairwise Plot
G2G_Res <- readRDS('../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/G2G_results.rds')
#Parse G2G Results into dataframe
G2G_Res <- G2G_Res$ALL
G2G_Res <- G2G_Res[sapply(G2G_Res,function(x) nrow(x) > 0)]

G2G_Res_Expand <- data.frame(Host_Chr = character(),Host_Pos = numeric(),Mtb_Pos = numeric(),Mtb_Gene=character(),P = numeric())
for(i in 1:length(G2G_Res)){
  cur_mtb_variant <- names(G2G_Res)[i]
  Host_Pos <- as.numeric(G2G_Res[[i]]$POS)
  Host_Chr  <- G2G_Res[[i]]$`#CHROM`
  Mtb_Pos <- as.numeric(strsplit(names(G2G_Res)[i],split = ':')[[1]][2])
  Mtb_Gene <- strsplit(names(G2G_Res)[i],split = ':')[[1]][1]
  P <- G2G_Res[[i]]$P

  G2G_Res_Expand <- rbind(G2G_Res_Expand,data.frame(Host_Chr=Host_Chr,Host_Pos=Host_Pos,Mtb_Pos=Mtb_Pos,Mtb_Gene=Mtb_Gene,P = P))
}

#Host Chr Length
host_chr_length <- getChromInfoFromUCSC("hg19")
host_chr_length <- dplyr::filter(host_chr_length,chrom %in% paste0('chr',c(1:22,'X'))) %>%
  dplyr::select(chr=chrom,size) %>%
  dplyr::mutate(chr = gsub(chr,pattern = 'chr',replacement = ''))
#Target Length
target_chr_length <- data.frame(chr='H37Rv',size = 4411532)

scale_factor <- sum(host_chr_length$size)/target_chr_length$size
target_chr_length <- data.frame(chr='H37Rv',size = 4411532*scale_factor)

#Change to Block Sturcture (For draw pairwise)
G2G_Block <- G2G_Res_Expand %>% dplyr::select(chr = Host_Chr,start = Host_Pos,tarSt=Mtb_Pos,P) %>%
  dplyr::mutate(orient=1,tarChr='H37Rv',tarSpecies = 'Mtb',tarSt=tarSt*scale_factor)


draw.pairwise <- function(dataTMP,ref_sizes,tar_sizes,refName,tarName,p_thresh) {
  xstart<-refchr<-tarchr<-x<-y<-group<-fill<-NULL
  data <- dplyr::select(dataTMP,tarchr=tarChr,tarstart=tarSt,refchr=chr,refstart=start,P)
  colnames(ref_sizes) = c("refchr", "size")
  colnames(tar_sizes) = c("tarchr", "size")

  #This adds gap in between reference chromosomes and convert to "linear" genome
  for (i in c(1:nrow(ref_sizes))){
    #print(i)
    if (i == 1){
      total_start = 1
      total_end = ref_sizes[i, "size"]
    } else {
      total_start = total_end + 6000000
      total_end = total_start + ref_sizes[i, "size"]
    }
    ref_sizes[i,"xstart"] = total_start
    ref_sizes[i, "xend"] = total_end
  }

  #This adds gap in between target chromosomes
  for (i in c(1:nrow(tar_sizes))){
    #print(i)
    if (i == 1){
      total_start = 1
      total_end = tar_sizes[i, "size"]
    } else {
      total_start = total_end + 6000000
      total_end = total_start + tar_sizes[i, "size"]
    }
    tar_sizes[i,"xstart"] = total_start
    tar_sizes[i, "xend"] = total_end
  }

  #This converts coordinates to linear genome and creates synteny polygon coordinates
  synteny = data.frame()
  for (i in c(1:nrow(data))){
    tar_chr = data[i,"tarchr"]
    ref_chr = data[i,"refchr"]
    dir = data[i, "dir"]

    tar_add = tar_sizes[as.character(tar_sizes$tarchr)==as.character(tar_chr),]$xstart
    ref_add = ref_sizes[as.character(ref_sizes$refchr)==as.character(ref_chr),]$xstart
    tar_y = 0.1
    ref_y = 2
    tar_xstart = data[i,"tarstart"] + tar_add
    ref_xstart = data[i,"refstart"] + ref_add

    df = data.frame(x = ref_xstart,y=ref_y,xend=tar_xstart,yend=tar_y,Sig = as.factor(ifelse(data$P[i] < p_thresh,'Sig','NonSig')))
    synteny = rbind(synteny,df)
  }

  #making sure chr columns are factors
  tar_sizes$tarchr<-as.factor(tar_sizes$tarchr)
  ref_sizes$refchr<-as.factor(ref_sizes$refchr)


  #This prints plot
  print(ggplot2::ggplot(size = 0.5, font = 10, data = data) +
          ggplot2::geom_rect(data=ref_sizes, mapping=ggplot2::aes(xmin=xstart, xmax=xend, ymin=2, ymax=2.10, fill=refchr),
                             color="black", alpha = 0.85, size = 0.2 ) +
          ggplot2::geom_text(data=ref_sizes,ggplot2::aes(x=(xstart+xend)/2,y=2.25,label=refchr),size=4,angle = 45) +
          ggplot2::geom_text(mapping=ggplot2::aes(x=10,y=2.4, label=refName),size=5) +
          ggplot2::geom_rect(data=tar_sizes, mapping=ggplot2::aes(xmin=xstart, xmax=xend, ymin=0, ymax=0.10),fill="purple",
                             color="black", alpha = 0.85, size = 0.2 ) +
          ggplot2::geom_text(data=tar_sizes,ggplot2::aes(x=(xstart+xend)/2,y=-0.2,label=tarchr),size=4) +
          ggplot2::geom_text(mapping=ggplot2::aes(x=10,y=-0.05, label=tarName),size=5) +
          ggplot2::geom_segment(data = synteny, ggplot2::aes(x = x, y = y,xend=xend,yend=yend,colour=Sig,alpha=Sig)) +
          ggplot2::scale_fill_manual(values = c("1" = "#BFD73B", "2" = "#39ACE2", "3" = "#F16E8A",
                                                "4" = "#2DB995", "5" = "#855823", "6" = "#A085BD",
                                                "7" = "#2EB560", "8" = "#D79128", "9" = "#FDBB63",
                                                "10" = "#AFDFE5", "11" = "#BF1E2D", "12" = "purple4",
                                                "13"= "#B59F31", "14" = "#F68B1F", "15" = "#EF374B",
                                                "16" = "#D376FF", "17" = "#009445", "18" = "#CE4699",
                                                "19" = "#7C9ACD", "20" = "#84C441", "21" = "#404F23",
                                                "22" = "#607F4B",  "X" = "blue")) +
          ggplot2::scale_color_manual(values=c('NonSig'='grey85','Sig'='red'))+
          ggplot2::scale_alpha_manual(values=c('NonSig'=0.3,'Sig'=1)) +
          ggplot2::theme(panel.background = ggplot2::element_blank(),
                         strip.background = ggplot2::element_blank(),
                         axis.title.y = ggplot2::element_blank(),
                         axis.title.x = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_blank(),
                         axis.text.y = ggplot2::element_blank(),
                         axis.ticks.x=ggplot2::element_blank(),
                         axis.ticks.y=ggplot2::element_blank(),
                         legend.position="none")
  )


}
P_Thresh <- 1.02459e-10
draw.pairwise(rbind(G2G_Block  %>% dplyr::filter(P >= P_Thresh),G2G_Block %>% dplyr::filter(P < P_Thresh)),host_chr_length,target_chr_length,'Human','M.tb',1.02459e-10)

#Pathogen Plot
Sig_Hits <- lapply(G2G_Res,function(x) dplyr::filter(x,P <= P_Thresh))
Sig_Hits <- Sig_Hits[sapply(Sig_Hits,function(x) nrow(x) > 0)]
Sig_Hits_Pos <- as.numeric(sapply(names(Sig_Hits),function(x) strsplit(x=x,split = ':')[[1]][2]))

results_dir <- '../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/'
Mtb_Variants <- dir(results_dir)
Mtb_Variants <- Mtb_Variants[grepl(Mtb_Variants,pattern = 'hybrid.gz')]
Mtb_Variants_Pos <- as.numeric(sapply(Mtb_Variants,function(x) strsplit(x=x,split = ':')[[1]][2]))

library(gggenes)
H37Rv_GFF <- data.table::fread('../../data/Mtb/H37Rv.gff3',skip = 2,sep = '\t') %>% dplyr::select(Type=V3,Start = V4,End=V5,orientation=V7,Annot=V9) %>% dplyr::filter(Type=='gene') %>%
  dplyr::select(-Type) %>% dplyr::mutate(forward = ifelse(orientation == '+',T,F))
H37Rv_GFF$Gene_Name <- sapply(H37Rv_GFF$Annot,function(x) strsplit(strsplit(x=x,split = 'Name=')[[1]][2],split = ';')[[1]][1])
H37Rv_GFF$Locus_Tag <- sapply(H37Rv_GFF$Annot,function(x) strsplit(strsplit(x=x,split = 'locus_tag=')[[1]][2],split = ';')[[1]][1])
H37Rv_GFF <- H37Rv_GFF %>% dplyr::select(-Annot)

Pos_Offset <- 5000
pathogen_plots <- list()
for(i in 1:length(Sig_Hits)){
  Neighbours <- Mtb_Variants_Pos[(Mtb_Variants_Pos > Sig_Hits_Pos[i] - Pos_Offset) & (Mtb_Variants_Pos < Sig_Hits_Pos[i] + Pos_Offset)]
  Neighbours_df <- lapply(Mtb_Variants[match(Neighbours,Mtb_Variants_Pos)],function(x) data.table::fread(cmd = glue::glue('zcat {results_dir}{x}')))
  Neighbours_df <- lapply(Neighbours_df,function(x) dplyr::filter(x,ID == Sig_Hits[[i]]$ID))
  Neighbours_df <- data.frame(POS = Neighbours,P=sapply(Neighbours_df,function(x) x$P),ID = sapply(Mtb_Variants[match(Neighbours,Mtb_Variants_Pos)],function(x) paste0(strsplit(paste0(strsplit(x,split = ':')[[1]][c(1,3)],collapse = '_'),split='\\.')[[1]][1:2],collapse = '.')))
  Neighbours_Genes <- H37Rv_GFF %>% dplyr::filter(Start > Sig_Hits_Pos[i] - Pos_Offset,End < Sig_Hits_Pos[i] + Pos_Offset)
  pathogen_plots[[i]] <- ggplot() +
    geom_gene_arrow(aes(xmin = Start, xmax = End,y=-1.5, fill = Gene_Name,label = Gene_Name,forward = forward),arrowhead_height = unit(10, "mm"), arrowhead_width = unit(3, "mm"),arrow_body_height = unit(10,'mm'),data = Neighbours_Genes) +
    geom_gene_label(aes(xmin = Start, xmax = End,y=-1.5, fill = Gene_Name,label = Gene_Name,forward = forward),align = "left",data = Neighbours_Genes) + geom_point(aes(x=POS,y=-log10(P)),data = Neighbours_df) + ylim(-2,12) +
    labs(y=TeX("$-log_{10}(P)$"),x= 'M.tb Genome Position') + ggtitle(Sig_Hits[[i]]$ID) +
    geom_text_repel(aes(x=POS,y=-log10(P),label = ID),data = dplyr::filter(Neighbours_df,P == min(Neighbours_df$P)))

}

#Host Plot
hg19_GFF <- data.table::fread(cmd = "cat ../../data/Genotyping/gencode.v19.annotation.gff3 | grep -v '#'",sep = '\t') %>% dplyr::select(Chr=V1,Type=V3,Start = V4,End=V5,orientation=V7,Annot=V9) %>% dplyr::filter(Type=='gene') %>%
  dplyr::select(-Type) %>% dplyr::mutate(forward = ifelse(orientation == '+',T,F))
hg19_GFF$Gene_Name <- sapply(hg19_GFF$Annot,function(x) strsplit(strsplit(x=x,split = 'gene_name=')[[1]][2],split = ';')[[1]][1])
hg19_GFF$Gene_Type <- sapply(hg19_GFF$Annot,function(x) strsplit(strsplit(x=x,split = 'gene_type=')[[1]][2],split = ';')[[1]][1])
hg19_GFF <- hg19_GFF %>% dplyr::filter(Gene_Type == 'protein_coding') %>% dplyr::select(-Annot,-Gene_Type)

Pos_Offset <- 200000
host_plots <- list()
for(i in 1:length(Sig_Hits)){
  GWAS_Data <- data.table::fread(cmd = glue::glue('zcat {results_dir}{names(Sig_Hits)[i]}'))
  GWAS_Data <- GWAS_Data %>% dplyr::filter(`#CHROM` == Sig_Hits[[i]]$`#CHROM`) %>%
    dplyr::filter(POS > Sig_Hits[[i]]$POS - Pos_Offset,POS < Sig_Hits[[i]]$POS + Pos_Offset) %>% dplyr::select(POS,ID,P)
  #Get LD Matrix
  system(glue::glue('~/Software/plink --bfile ../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/LINEAGE_ALL/TB_DAR_G2G --r2 inter-chr --ld-snp {GWAS_Data$ID[which.min(GWAS_Data$P)]} --ld-window-r2 0 --out ../../scratch/tmp'))
  ld_mat <- data.table::fread('../../scratch/tmp.ld') %>% dplyr::filter(SNP_B %in% GWAS_Data$ID) %>%
    dplyr::select(ID=SNP_B,r2=R2)
  GWAS_Data <- dplyr::left_join(GWAS_Data,ld_mat,by=c('ID'='ID'))
  system('rm ../../scratch/tmp.ld')
  Neighbours_Genes <- hg19_GFF %>% dplyr::filter(Chr == paste0('chr',Sig_Hits[[i]]$`#CHROM`),Start > Sig_Hits[[i]]$POS - Pos_Offset,End < Sig_Hits[[i]]$POS + Pos_Offset)
  host_plots[[i]] <- ggplot() +
    geom_gene_arrow(aes(xmin = Start, xmax = End,y=ifelse(forward,-1.5,-2.5), fill = Gene_Name,label = Gene_Name,forward = forward),arrowhead_height = unit(10, "mm"), arrowhead_width = unit(3, "mm"),arrow_body_height = unit(10,'mm'),data = Neighbours_Genes) +
    geom_gene_label(aes(xmin = Start, xmax = End,y=ifelse(forward,-1.5,-2.5), fill = Gene_Name,label = Gene_Name,forward = forward),align = "left",data = Neighbours_Genes) + geom_point(aes(x=POS,y=-log10(P),color=r2),data = GWAS_Data) + ylim(-3,12) +
    labs(y=TeX("$-log_{10}(P)$"),x= paste0('Human Chr',Sig_Hits[[i]]$`#CHROM`,' Position')) + ggtitle(paste0(strsplit(paste0(strsplit(names(Sig_Hits)[i],split = ':')[[1]][c(1,3)],collapse = '_'),split='\\.')[[1]][1:2],collapse = '.')) +
    geom_text_repel(aes(x=POS,y=-log10(P),label = ID),data = dplyr::filter(GWAS_Data,P == min(GWAS_Data$P)))

}

#QQ Plot
Sig_Hits_Sum_Stats <- lapply(names(Sig_Hits),function(x) data.table::fread(glue::glue("zcat {results_dir}{names(Sig_Hits)[i]}")))
qqman::qq(Sig_Hits_Sum_Stats[[1]]$P)
text(paste0('Lambda=',round(GWAS.utils::genomic_inflation(P=Sig_Hits_Sum_Stats[[1]]$P),2)),x=2,y = 8)
qqman::qq(Sig_Hits_Sum_Stats[[2]]$P)
text(paste0('Lambda=',round(GWAS.utils::genomic_inflation(P=Sig_Hits_Sum_Stats[[2]]$P),2)),x=2,y = 8)

#Summary Statistics (G2G Hits)
bac_load <- data.table::fread('../data/pheno/metadata_Sinergia_final_dataset_human_bac_genome_available_corrected.txt') %>% dplyr::select(PATIENT_ID,Bacterial_load)
G2G_Obj <- readRDS('../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds')
AA_Matrix <- as.data.frame(G2G_Obj$aa_matrix_full[,c("Rv2348c_Rv2348c:2626678:p.Ile101Met",'fixA_Rv3029c:3388671:p.Thr67Met')])
AA_Matrix[AA_Matrix==2] <- 1

G2G_Pheno <- cbind(cbind(cbind(cbind(G2G_Obj$covars$ALL %>% dplyr::select(patient_sex,HIV_status,TB_RF_smoking,age,BMI),AA_Matrix),
                   G2G_Obj$host_PCs$ALL %>% dplyr::select(PC1,PC2,PC3)),G2G_Obj$tb_score$ALL %>% dplyr::select(TB_score)),G2G_Obj$xray_score$ALL %>% dplyr::select(Xray_score))
G2G_Pheno$HIV_status[G2G_Pheno$HIV_status=='NA'] <- NA
G2G_Pheno$age <- scale(as.numeric(G2G_Pheno$age))
G2G_Pheno$BMI <- scale(as.numeric(G2G_Pheno$BMI))
G2G_Pheno$TB_RF_smoking <- ifelse(G2G_Pheno$TB_RF_smoking=='yes',1,0)

fixA_mdl_sex <- summary(glm(formula = '`fixA_Rv3029c:3388671:p.Thr67Met`~patient_sex',data = G2G_Pheno))$coefficients['patient_sexmale',c('Estimate','Std. Error','Pr(>|t|)')]
fixA_mdl_age <- summary(glm(formula = '`fixA_Rv3029c:3388671:p.Thr67Met`~ age',data = G2G_Pheno))$coefficients['age',c('Estimate','Std. Error','Pr(>|t|)')]
fixA_mdl_PC <- summary(glm(formula = '`fixA_Rv3029c:3388671:p.Thr67Met`~PC1+PC2+PC3',data = G2G_Pheno))$coefficients[c('PC1','PC2','PC3'),c('Estimate','Std. Error','Pr(>|t|)')]
fixA_mdl_HIV <- summary(glm(formula = '`fixA_Rv3029c:3388671:p.Thr67Met`~HIV_status',data = G2G_Pheno))$coefficients['HIV_statusnegative',c('Estimate','Std. Error','Pr(>|t|)')]
fixA_mdl_Smoking <- summary(glm(formula = '`fixA_Rv3029c:3388671:p.Thr67Met`~TB_RF_smoking',data = G2G_Pheno))$coefficients['TB_RF_smoking',c('Estimate','Std. Error','Pr(>|t|)')]
fixA_mdl_BMI <- summary(glm(formula = '`fixA_Rv3029c:3388671:p.Thr67Met`~BMI',data = G2G_Pheno))$coefficients['BMI',c('Estimate','Std. Error','Pr(>|t|)')]

Rv2348c_mdl_sex <- summary(glm(formula = '`Rv2348c_Rv2348c:2626678:p.Ile101Met`~patient_sex',data = G2G_Pheno))$coefficients['patient_sexmale',c('Estimate','Std. Error','Pr(>|t|)')]
Rv2348c_mdl_age <- summary(glm(formula = '`Rv2348c_Rv2348c:2626678:p.Ile101Met`~ age',data = G2G_Pheno))$coefficients['age',c('Estimate','Std. Error','Pr(>|t|)')]
Rv2348c_mdl_PC <- summary(glm(formula = '`Rv2348c_Rv2348c:2626678:p.Ile101Met`~PC1+PC2+PC3',data = G2G_Pheno))$coefficients[c('PC1','PC2','PC3'),c('Estimate','Std. Error','Pr(>|t|)')]
Rv2348c_mdl_HIV <- summary(glm(formula = '`Rv2348c_Rv2348c:2626678:p.Ile101Met`~HIV_status',data = G2G_Pheno))$coefficients['HIV_statusnegative',c('Estimate','Std. Error','Pr(>|t|)')]
Rv2348c_mdl_Smoking <- summary(glm(formula = '`Rv2348c_Rv2348c:2626678:p.Ile101Met`~TB_RF_smoking',data = G2G_Pheno))$coefficients['TB_RF_smoking',c('Estimate','Std. Error','Pr(>|t|)')]
Rv2348c_mdl_BMI <- summary(glm(formula = '`Rv2348c_Rv2348c:2626678:p.Ile101Met`~BMI',data = G2G_Pheno))$coefficients['BMI',c('Estimate','Std. Error','Pr(>|t|)')]

#Patient Characteristics (Lineage)
library(table1)
pheno <- data.table::fread('../../data/pheno/metadata_Sinergia_final_dataset_human_bac_genome_available_corrected.txt')
pheno$HIV_status[pheno$HIV_status==""] <- NA


host_PC <- data.table::fread('../../scratch/Host/TB_DAR_GWAS_PCA.eigenvec') %>%
  dplyr::select(IID=V1,PC1=V3,PC2=V4,PC3=V5,PC4=V6)
G2G_Obj <- readRDS('../../scratch/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/G2G_Obj.rds')


host_PC <- host_PC %>% dplyr::inner_join(G2G_Obj$both_IDs_to_keep$ALL,by=c('IID'='FAM_ID')) %>% dplyr::select(PATIENT_ID,PC1,PC2,PC3,PC4)

pheno <- pheno %>% dplyr::left_join(host_PC)

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.0f %%)", FREQ, PCT))))
}
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- round(summary(aov(y ~ g))[[1]][["Pr(>F)"]][1],2)
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- round(chisq.test(table(y, g))$p.value,2)
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

table1(~ factor(patient_sex) + age + BMI + TB_RF_smoking + HIV_status + TB_score + Ct_value + Xray_score | Lineage, 
       data=pheno,render.continuous=my.render.cont, 
       render.categorical=my.render.cat,extra.col=list(`P-value`=pvalue))

host_PC_eigenval <- data.table::fread('../../scratch/Host/TB_DAR_GWAS_PCA.eigenval')
ggplot2::ggplot(pheno) + geom_point(aes(x=PC1,y=PC2,color = Lineage)) + 
  labs(x = paste0('PC1 (',round(host_PC_eigenval$V1[1] / sum(host_PC_eigenval$V1) * 100,2),'%)'),y=paste0('PC2 (',round(host_PC_eigenval$V1[2] / sum(host_PC_eigenval$V1) * 100,2),'%)'))
ggplot2::ggplot(pheno) + geom_point(aes(x=PC3,y=PC4,color = Lineage)) + 
  labs(x = paste0('PC3 (',round(host_PC_eigenval$V1[3] / sum(host_PC_eigenval$V1) * 100,2),'%)'),y=paste0('PC4 (',round(host_PC_eigenval$V1[4] / sum(host_PC_eigenval$V1) * 100,2),'%)'))



