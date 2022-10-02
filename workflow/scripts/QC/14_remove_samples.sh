#!/bin/bash
QC_DIR="../../../data/Genotyping_WGS/QC/"
cat $QC_DIR"Genetic_Duplicates.txt" | grep -v 'IID' > $QC_DIR"tmp1.txt"
cat $QC_DIR"PCA_Outliers.txt" | grep -v 'IID' >> $QC_DIR"tmp1.txt"
cat $QC_DIR"Regenotyped_Excl.txt" | grep -v 'IID' >> $QC_DIR"tmp1.txt"
cat $QC_DIR"Relatives.txt" | grep -v 'IID' >> $QC_DIR"tmp1.txt"
cat $QC_DIR"Sex_Mismatch.txt" | grep -v 'IID' >> $QC_DIR"tmp1.txt"
echo "FID IID" > $QC_DIR"All_Excl_Samples.txt"
sort $QC_DIR"tmp1.txt" | uniq >> $QC_DIR"All_Excl_Samples.txt"
rm $QC_DIR"tmp1.txt"
~/Software/plink2 --vcf ../../../data/Genotyping_WGS/TBDAR.WGS.Imputed.AnalysisReady.vcf.gz --remove $QC_DIR"All_Excl_Samples.txt" --make-bed  --double-id --out ../../../data/Genotyping_WGS/TBDAR.WGS.Imputed.GWASReady
