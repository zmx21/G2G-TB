source('CalcM_eff.R')

Variant_Set_Table <- data.frame(Mtb_Variants = sapply(Variant_Sets,function(x) paste0(x,collapse = ';')),Variant_Set = paste0('Variant_Set_',1:length(Variant_Sets)))


out_dir <- '../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/Summary_Stats/'
system(glue::glue("mkdir -p {out_dir}"))
system(glue::glue("mkdir -p {out_dir}Data/"))

results_dir <- '../../results/Burden_False_SIFT_False_Del_False_HomoOnly_True_HetThresh_10/PLINK/PC_3_pPC_0/Stratified_False/LINEAGE_ALL/'

for(i in 1:length(Variant_Sets)){
  print(i)
  #take first variant of variant set
  cur_variant <- Variant_Sets[[i]][1]
  cur_variant_set <- Variant_Set_Table$Variant_Set[i]
  #path to full PLINK results
  cur_summary_stats <- glue::glue("{results_dir}{cur_variant}.{cur_variant}.glm.logistic.hybrid.gz")
  
  #Write out summary stats
  system(glue::glue("zcat {cur_summary_stats} | awk '{{print $3,$10,$11,$13}}' > {out_dir}Data/{cur_variant_set}.txt"))
  system(glue::glue("pigz {out_dir}Data/{cur_variant_set}.txt"))
  
}

#write pathogen variant table
data.table::fwrite(Variant_Set_Table,file=glue::glue("{out_dir}Mtb_Variants_Table.txt"),sep = ' ',quote = F)
#write host variant table
system(glue::glue("zcat {cur_summary_stats} | awk '{{print $1,$2,$3,$4,$5}}' > {out_dir}Human_Variants_Table.txt"))
