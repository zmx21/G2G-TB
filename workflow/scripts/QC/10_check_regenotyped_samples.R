TB_DAR_Samples <- system('~/Software/bcftools query -l ../../../Genotyping_WGS/TBDAR.WGS.Imputed.AnalysisReady.vcf.gz',intern = T)
Patient_IDs <- sapply(TB_DAR_Samples,function(x) strsplit(x,split = '_')[[1]][2])
Patient_IDs <- sapply(Patient_IDs,function(x) gsub(x=x,pattern = 'Fellay.',replacement = ''))

Dup_IDs <- TB_DAR_Samples[which(Patient_IDs %in% Patient_IDs[duplicated(Patient_IDs)])]
Dup_IDs_df <- data.frame(FID = Dup_IDs,IID = Dup_IDs)
data.table::fwrite(Dup_IDs_df,'../../../Genotyping_WGS/regenotyped.ids',sep = ' ',row.names = F,quote = F)
system('~/Software/plink2 --bfile ../../../Genotyping_WGS/TBDAR.WGS.Imputed.AnalysisReady --chr 1-22 --keep ../../../Genotyping_WGS/regenotyped.ids --make-bed --out ../../../Genotyping_WGS/TBDAR.Regenotyped')

regenotyped_plink <- snpStats::read.plink('../../../Genotyping_WGS/TBDAR.Regenotyped')

WGS_Genotypes <- regenotyped_plink$genotypes[sapply(regenotyped_plink$fam$member,function(x) grepl(x=x,pattern = 'WGS')),]
rownames(WGS_Genotypes) <- sapply(rownames(WGS_Genotypes) ,function(x) strsplit(x=x,split = '\\.')[[1]][2])

H3Africa_Genotypes <- regenotyped_plink$genotypes[!sapply(regenotyped_plink$fam$member,function(x) grepl(x=x,pattern = 'WGS')),]
H3Africa_Genotypes <- H3Africa_Genotypes[!duplicated(sapply(rownames(H3Africa_Genotypes),function(x) strsplit(x=x,split = '_')[[1]][2])),]
rownames(H3Africa_Genotypes) <- sapply(rownames(H3Africa_Genotypes),function(x) strsplit(x=x,split = '_')[[1]][2])
H3Africa_Genotypes <- H3Africa_Genotypes[rownames(WGS_Genotypes),]

snp_comp <- snpStats::sm.compare(H3Africa_Genotypes,WGS_Genotypes,col.wise = T)

sample_diff <- (snp_comp$row.wise[,'Disagree'] - snp_comp$row.wise[,'NA.disagree']) / (snp_comp$row.wise[,'Agree'] + snp_comp$row.wise[,'Disagree'] - snp_comp$row.wise[,'NA.disagree'])
snp_diff <- (snp_comp$col.wise[,'Disagree'] - snp_comp$col.wise[,'NA.disagree']) / (snp_comp$col.wise[,'Agree'] + snp_comp$col.wise[,'Disagree'] - snp_comp$col.wise[,'NA.disagree'])
hist(sample_diff,main = 'Percent of SNPs Different per Sample',xlab = 'Percent')
hist(snp_diff,main = 'Percent of Samples Different per SNP',xlab = 'Percent')

#Exclude duplicate between batch 1 and batch2. Exclude WGS Samples
WGS_Dup_IDs <- dplyr::filter(Dup_IDs_df,grepl(x = IID,pattern = 'WGS'))
H3Africa_Dup_IDs <- dplyr::filter(Dup_IDs_df,!grepl(x = IID,pattern = 'WGS'))
H3Africa_Dup_IDs$Simple_ID <- sapply(H3Africa_Dup_IDs$IID,function(x) strsplit(x=x,split = '_')[[1]][2])

#Find mismatch Duplicates
mismatching_Dup_IDs <- regenotyped_plink$fam$member[sapply(regenotyped_plink$fam$member,function(x) any(sapply(names(sample_diff[sample_diff > 0.1]),function(y) grepl(x=x,pattern = y))))]
data.table::fwrite(data.frame(FID = mismatching_Dup_IDs,IID=mismatching_Dup_IDs),'../../../Genotyping_WGS/QC/Regenotyped_Discordant.txt',sep = ' ',quote = F,row.names = F)

Dup_IDs_to_excl <- rbind(rbind(WGS_Dup_IDs,H3Africa_Dup_IDs[duplicated(H3Africa_Dup_IDs$Simple_ID),] %>% dplyr::select(FID,IID)),data.frame(IID= mismatching_Dup_IDs,FID=mismatching_Dup_IDs)) %>% dplyr::distinct(FID,IID,.keep_all = T)
data.table::fwrite(Dup_IDs_to_excl,'../../../Genotyping_WGS/QC/Regenotyped_Excl.txt',row.names = F,sep = ' ',quote = F)
