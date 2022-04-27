#!/bin/bash
#Input Files
IMPUTED_VCF='../Imputed/AFGR_Tanz/TBDAR.Imputed.AFGR.Tanz.vcf.gz'
INFO_FILE="../Imputed/AFGR_Tanz/WGS_AFGR_Merged_Info.txt"
WGS_VCF='../../../WGS/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz'
#Thresholds
MAF=0.05
HWE=1e-5
GENO=0.1
MIND=0.1
INFO_SCORE=0.8
#Parameters
N_CORES=50
OUT_DIR='../../../Genotyping_WGS/'
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR'tmp/'

#Filter WGS VCF 
#~/Software/plink --vcf $WGS_VCF --threads $N_CORES --keep-allele-order --hwe $HWE --geno $GENO --mind $MIND --export vcf bgz --out $OUT_DIR"tmp/WGS.tmp1"
~/Software/bcftools index -t --threads $N_CORES $OUT_DIR"tmp/WGS.tmp1.vcf.gz"

#Filter Imputed VCF
#awk -v info=$INFO_SCORE '{if($6 >= info && NR != 1) print $3}' $INFO_FILE > $OUT_DIR"tmp/Imputed.snps.pass.info.txt"
#~/Software/plink --vcf $IMPUTED_VCF --extract $OUT_DIR"tmp/Imputed.snps.pass.info.txt" --threads $N_CORES --keep-allele-order --hwe $HWE --geno $GENO --mind $MIND --export vcf bgz --out $OUT_DIR"tmp/Imputed.tmp1"
#~/Software/bcftools index -t --threads $N_CORES $OUT_DIR"tmp/Imputed.tmp1.vcf.gz"
#Merge Imputed and WGS
#~/Software/bcftools merge -m none --threads $N_CORES -O z -o $OUT_DIR"TBDAR.WGS.Imputed.Raw.vcf.gz" $OUT_DIR"tmp/Imputed.tmp1.vcf.gz" $OUT_DIR"tmp/WGS.tmp1.vcf.gz"

#Apply Filter to Merged File
~/Software/plink --vcf $OUT_DIR"TBDAR.WGS.Imputed.Raw.vcf.gz" --threads $N_CORES --keep-allele-order --hwe $HWE --geno $GENO --maf $MAF --export vcf bgz --out $OUT_DIR"TBDAR.WGS.Imputed.AnalysisReady.vcf.gz"
