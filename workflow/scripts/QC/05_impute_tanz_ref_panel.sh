#!/bin/bash
CHR=$1
SOFTWARE_DIR=~/Software/
REF_DIR='../../../WGS/Ref_Panel/'
INPUT='../Filt_Merged/TBDAR.Genotyped.Merged.Filt'
OUT_DIR='../Imputed/Tanz/'
N_THREAD=5
mkdir -p $OUT_DIR"tmp"
if [ $CHR == 'X' ]
then
	echo $CHR
	#Process X CHR Autosomal Region
	for PAR in PAR1 PAR2 
	do

		#Extract PAR
		$SOFTWARE_DIR"plink2" --bfile $INPUT --export vcf --chr $PAR --out $OUT_DIR"tmp/chrX"$PAR
		$SOFTWARE_DIR"bcftools" norm -d none -O v -o $OUT_DIR"tmp/chrX"$PAR.nodup.vcf $OUT_DIR"tmp/chrX"$PAR.vcf
		#Phase
		$SOFTWARE_DIR"shapeit" -V $OUT_DIR"tmp/chrX"$PAR".nodup.vcf" \
			-M $REF_DIR"1000GP_Phase3/genetic_map_chrX_"$PAR"_combined_b37.txt" \
                	-O $OUT_DIR"tmp/chr"$PAR".phased" \
                	--thread $N_THREAD \
                	--force
        	$SOFTWARE_DIR"shapeit" -convert \
                	--input-haps $OUT_DIR"tmp/chr"$PAR".phased" \
                	--output-vcf $OUT_DIR"tmp/chr"$PAR".phased.vcf"
		bgzip $OUT_DIR"tmp/chr"$PAR".phased.vcf"
		$SOFTWARE_DIR"bcftools" index -t $OUT_DIR"tmp/chr"$PAR".phased.vcf.gz"
	done

	#Merge phased PAR
	$SOFTWARE_DIR"bcftools" concat -O z -o $OUT_DIR"tmp/chrXPAR.phased.vcf.gz" $OUT_DIR"tmp/chrPAR1.phased.vcf.gz" $OUT_DIR"tmp/chrPAR2.phased.vcf.gz"
	$SOFTWARE_DIR"bcftools" index -t $OUT_DIR"tmp/chrXPAR.phased.vcf.gz"

	#Impute
	$SOFTWARE_DIR"Minimac3-omp" --refHaps $REF_DIR"m3vcfs/chrXPAR.m3vcf.gz" \
		--haps $OUT_DIR"tmp/chrXPAR.phased.vcf.gz" \
		--prefix $OUT_DIR"chrXPAR"
	$SOTWARE_DIR"bcftools" index -t $OUT_DIR"chrXPAR.dose.vcf.gz"	
	
	#Process X CHR
	$SOFTWARE_DIR"plink2" --bfile $INPUT --export vcf --chr "X" --out $OUT_DIR"tmp/chrX"
	$SOFTWARE_DIR"bcftools" norm -d none -O v -o $OUT_DIR"tmp/chrX.nodup.vcf" $OUT_DIR"tmp/chrX.vcf"
	#Phase
	$SOFTWARE_DIR"shapeit" -V $OUT_DIR"tmp/chrX.nodup.vcf" \
		-M $REF_DIR"1000GP_Phase3/genetic_map_chrX_nonPAR_combined_b37.txt" \
                -O $OUT_DIR"tmp/chrX.phased" \
		--thread $N_THREAD \
                --force

        $SOFTWARE_DIR"shapeit" -convert \
		--input-haps $OUT_DIR"tmp/chrX.phased" \
		--output-vcf $OUT_DIR"tmp/chrX.phased.vcf"
        bgzip $OUT_DIR"tmp/chrX.phased.vcf"

	$SOFTWARE_DIR"minimac4" --refHaps $REF_DIR"m3vcfs/chrX.m3vcf.gz" \
		--haps $OUT_DIR"tmp/chrX.phased.vcf.gz" \
		--prefix $OUT_DIR"chrX" --chunkLengthMb 100 --minRatio 0.05
	$SOFTWARE_DIR"bcftools" index -t $OUT_DIR"chrX.dose.vcf.gz"		
else
	echo $CHR
	#Phase
	$SOFTWARE_DIR"plink2" --bfile $INPUT --export vcf --chr $CHR --out $OUT_DIR"tmp/chr"$CHR
	$SOFTWARE_DIR"bcftools" norm -d none -O v -o  $OUT_DIR"tmp/chr"$CHR".nodup.vcf"  $OUT_DIR"tmp/chr"$CHR".vcf"

	$SOFTWARE_DIR"shapeit" -V $OUT_DIR"tmp/chr"$CHR".nodup.vcf" \
		-M $REF_DIR"1000GP_Phase3/genetic_map_chr"$CHR"_combined_b37.txt" \
		-O $OUT_DIR"tmp/chr"$CHR".phased" \
                --thread $N_THREAD \
                --force

        $SOFTWARE_DIR"shapeit" -convert \
                --input-haps $OUT_DIR"tmp/chr"$CHR".phased" \
                --output-vcf $OUT_DIR"tmp/chr"$CHR".phased.vcf"
        bgzip $OUT_DIR"tmp/chr"$CHR".phased.vcf"
	$SOTWARE_DIR"bcftools" index -t $OUT_DIR"tmp/chr"$CHR".phased.vcf.gz"
	
	#Impute
	$SOFTWARE_DIR"minimac4" --refHaps $REF_DIR"m3vcfs/chr"$CHR".m3vcf.gz" \
		--haps $OUT_DIR"tmp/chr"$CHR".phased.vcf.gz" \
		--prefix $OUT_DIR"chr"$CHR
	$SOFTWARE_DIR"bcftools" index -t $OUT_DIR"chr"$CHR".dose.vcf.gz"	
fi
