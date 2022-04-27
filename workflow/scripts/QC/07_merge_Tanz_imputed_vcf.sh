#!/bin/bash
DATA_DIR='../Imputed/Tanz/'
OUTPUT='TBDAR.Tanz.Imputed'
~/Software/bcftools concat $DATA_DIR"chr"{1..22}".dose.vcf.gz" $DATA_DIR"chrX.dose.vcf.gz" -O z -o $DATA_DIR$OUTPUT".vcf.gz"
~/Software/bcftools index -t --threads 5 $DATA_DIR$OUTPUT".vcf.gz"
~/Software/bcftools query -f "%CHROM %POS %ID %REF %ALT %INFO/R2 %INFO/AF\n" $DATA_DIR$OUTPUT".vcf.gz" > $DATA_DIR$OUTPUT".info.txt"
~/Software/bcftools query -f "%CHROM %POS %ID %REF %ALT %INFO/R2 %INFO/AF\n" $DATA_DIR"chrXPAR.dose.vcf.gz" > $DATA_DIR"chrXPAR.info.txt"
