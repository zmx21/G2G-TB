#!/bin/bash
DATA_DIR='../Imputed/AFGR/'
OUTPUT='TBDAR.AFGR.Imputed'
#~/Software/bcftools concat $DATA_DIR""{1..22}".vcf.gz" $DATA_DIR"X.vcf.gz" -O z -o $DATA_DIR$OUTPUT".vcf.gz"
~/Software/bcftools index -t --threads 5 $DATA_DIR$OUTPUT".vcf.gz"
~/Software/bcftools query -f "%CHROM %POS %ID %REF %ALT %INFO/RefPanelAF %INFO/INFO\n" $DATA_DIR$OUTPUT".vcf.gz" > $DATA_DIR$OUTPUT".info.txt"
~/Software/plink2 --vcf $DATA_DIR$OUTPUT".vcf.gz" --max-alleles 2 --set-missing-var-ids "@:#[b37]\$r,\$a" --const-fid --make-bed --out $DATA_DIR$OUTPUT
cut -f 2 $DATA_DIR$OUTPUT".bim" | sort | uniq -d > $DATA_DIR$OUTPUT".duprs"
~/Software/plink2 --bfile $DATA_DIR$OUTPUT --exclude $DATA_DIR$OUTPUT".duprs" --make-bed --out $DATA_DIR$OUTPUT".nodup"
~/Software/plink2 --bfile $DATA_DIR$OUTPUT".nodup" --export vcf --out $DATA_DIR$OUTPUT".nodup"
bgzip -c $DATA_DIR$OUTPUT".nodup.vcf" > $DATA_DIR$OUTPUT".nodup.vcf.gz"
~/Software/bcftools index -t --threads 5 $DATA_DIR$OUTPUT".nodup.vcf.gz"
~/Software/bcftools query -f "%CHROM %POS %ID %REF %ALT\n" $DATA_DIR$OUTPUT".nodup.vcf.gz" > $DATA_DIR$OUTPUT".nodup.vcf.gz.pos"
