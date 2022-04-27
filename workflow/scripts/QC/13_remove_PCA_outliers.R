GetPCAOutliers <- function(PCA_results,target_POP){
  pop_tbl <- data.table::fread('~/1KG_PCA/1KG_data/20130606_g1k.ped') %>% dplyr::select(ID=`Individual ID`,Population)
  pop_info <- data.table::fread('~/1KG_PCA/1KG_data/20131219.populations.tsv') %>% dplyr::select(Population = `Population Code`,SuperPopulation = `Super Population`)
  pop_tbl <- dplyr::left_join(pop_tbl,pop_info)
  
  genotyping_pca <- data.table::fread(PCA_results) %>% dplyr::select(ID=V2,PC1=V3,PC2=V4)
  genotyping_pca <- genotyping_pca %>% dplyr::left_join(pop_tbl %>% dplyr::select(ID,SuperPopulation))
  
  genotyping_pca$SuperPopulation[is.na(genotyping_pca$SuperPopulation)] <- 'TARGET'
  
  library(class)
  non_target_set <- genotyping_pca %>% dplyr::filter(SuperPopulation != 'TARGET' & SuperPopulation != 'AMR') %>% dplyr::select(-SuperPopulation,-ID) #Train KNN based on 1KG 
  non_target_set_class <- genotyping_pca %>% dplyr::filter(SuperPopulation != 'TARGET' & SuperPopulation != 'AMR') %>% dplyr::select(SuperPopulation) #Labels for the Train set
  non_target_set_class <- factor(non_target_set_class$SuperPopulation)
  
  train_sample <- sample(1:nrow(non_target_set),size = 0.7*nrow(non_target_set),replace = F)
  train_set <- non_target_set[train_sample,]
  train_set_class <- non_target_set_class[train_sample]
  test_set <- non_target_set[-train_sample,]
  
  target_set <- genotyping_pca %>% dplyr::filter(SuperPopulation == 'TARGET') 
  target_set_pred <- target_set  %>% dplyr::select(-ID,-SuperPopulation)
  knn_test <- knn(train = train_set,test = test_set,cl=train_set_class,k=4)
  print(glue::glue('Testing Accuracy: {sum(knn_test == non_target_set_class[-train_sample]) / length(knn_test)}'))
  
  knn_target <- knn(train = non_target_set,test = target_set_pred,cl=non_target_set_class,k=4) 
  outliers <- target_set$ID[knn_target != target_POP]
  p1 <- ggplot2::ggplot(genotyping_pca %>% dplyr::filter(SuperPopulation != 'AMR')) + aes(x=PC1,y=PC2,color = SuperPopulation) + geom_point() + ggrepel::geom_text_repel(data = genotyping_pca %>% dplyr::filter(ID %in%outliers),aes(x=PC1,y=PC2,label = ID))
  
  return(list(outliers=outliers,p1=p1))
  
}
# genotyping_outliers <- GetPCAOutliers('../../PCA/TBDAR.pca.eigenvec','AFR')
# WGS_outliers <- GetPCAOutliers('../../../WGS/PCA_WGS/TBDAR.WGS.eigenvec','AFR')

Genotyping_WGS_outliers <- GetPCAOutliers('../../../Genotyping_WGS/PCA/TBDAR.WGS.Imputed.eigenvec','AFR')
data.table::fwrite(data.frame(FID=Genotyping_WGS_outliers$outliers,IID=Genotyping_WGS_outliers$outliers),'../../../Genotyping_WGS/QC/PCA_Outliers.txt',sep = ' ')
