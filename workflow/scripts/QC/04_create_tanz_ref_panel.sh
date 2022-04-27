#!/bin/bash
INPUT_PREFIX='joined.hg19.nodup.nomismap'
REF_DIR="../../../WGS/Ref_Panel/1000GP_Phase3/"
DATA_DIR="../../../WGS/WGS_Host_Data/"
PHASED_DIR='../../../WGS/Ref_Panel/phased_vcfs/'
OUTPUT_DIR='../../../WGS/Ref_Panel/'

SOFTWARE_DIR=~/Software/
CHR=$1
NO_MISSING_OUT=$OUTPUT_DIR$INPUT_PREFIX".nomissing"

Geno_Thresh=0.9
N_THREAD=23

mkdir -p $PHASED_DIR
mkdir -p $OUTPUT_DIR"m3vcfs"

if [ $CHR == "X" ]; then
	#Run Check Sex
	#$SOFTWARE_DIR"plink" --bfile $DATA_DIR$INPUT_PREFIX --check-sex --out $NO_MISSING_OUT
	#awk '{ if($5=="PROBLEM") {print $1" "$2}}' $NO_MISSING_OUT".sexcheck" > $OUT_DIR"sex_mismatch.txt"
        
	#Process PAR like autosomes
	for PAR in PAR1 PAR2 
	do
		#Apply QC 
        	$SOFTWARE_DIR"plink2" --bfile $DATA_DIR$INPUT_PREFIX --split-par 'hg19' --geno $Geno_Thresh --chr $PAR --export vcf --out $NO_MISSING_OUT".chr"$PAR
        	#Run Phasing
        	$SOFTWARE_DIR"shapeit" -V $NO_MISSING_OUT".chr"$PAR".vcf" \
                	-M $REF_DIR"genetic_map_chrX_"$PAR"_combined_b37.txt" \
                	-O $NO_MISSING_OUT".chr"$PAR".phased" \
                	--thread $N_THREAD \
                	--force
        	$SOFTWARE_DIR"shapeit" -convert \
                	--input-haps $NO_MISSING_OUT".chr"$PAR".phased" \
                	--output-vcf $PHASED_DIR$INPUT_PREFIX".nomissing.chr"$PAR".phased.vcf"
		bgzip $PHASED_DIR$INPUT_PREFIX".nomissing.chr"$PAR".phased.vcf"
	done
	
	#Process X Chr
	#Apply QC
	$SOFTWARE_DIR"plink2" --bfile $DATA_DIR$INPUT_PREFIX --split-par 'hg19' --geno $Geno_Thresh --chr $CHR --export vcf --out $NO_MISSING_OUT".chr"$CHR

	#Run Phasing
	$SOFTWARE_DIR"shapeit" -V $NO_MISSING_OUT".chr"$CHR".vcf" \
		-M $REF_DIR"genetic_map_chrX_nonPAR_combined_b37.txt" \
		-O $NO_MISSING_OUT".chr"$CHR".phased" \
		--chrX \
		--thread $N_THREAD \
		--force
	$SOFTWARE_DIR"shapeit" -convert \
		--input-haps $NO_MISSING_OUT".chr"$CHR".phased" \
		--output-vcf $PHASED_DIR$INPUT_PREFIX".nomissing.chr"$CHR".phased.vcf"
	bgzip $PHASED_DIR$INPUT_PREFIX".nomissing.chr"$CHR".phased.vcf"

	#Create m3vcf for chrX
	$SOFTWARE_DIR"Minimac3-omp" --cpus 20 --refHap $PHASED_DIR$INPUT_PREFIX".nomissing.chrX.phased.vcf.gz" \
		--processReference --prefix $OUTPUT_DIR"m3vcfs/chrX"
	#Create m3vcf for PAR of X
	$SOFTWARE_DIR"bcftools" index -t $PHASED_DIR$INPUT_PREFIX".nomissing.chrPAR1.phased.vcf.gz"
	$SOFTWARE_DIR"bcftools" index -t $PHASED_DIR$INPUT_PREFIX".nomissing.chrPAR2.phased.vcf.gz"
	$SOFTWARE_DIR"bcftools" concat -O z -o $PHASED_DIR$INPUT_PREFIX".nomissing.chrXPAR.phased.vcf.gz" \
		$PHASED_DIR$INPUT_PREFIX".nomissing.chrPAR1.phased.vcf.gz" \
		$PHASED_DIR$INPUT_PREFIX".nomissing.chrPAR2.phased.vcf.gz" 
	$SOFTWARE_DIR"bcftools" index -t --threads 5 $PHASED_DIR$INPUT_PREFIX".nomissing.chrXPAR.phased.vcf.gz" 
	$SOFTWARE_DIR"Minimac3-omp" --cpus 20 --refHap $PHASED_DIR$INPUT_PREFIX".nomissing.chrXPAR.phased.vcf.gz" \
		--processReference --prefix $OUTPUT_DIR"m3vcfs/chrXPAR"

else
	#Apply QC
	$SOFTWARE_DIR"plink2" --bfile $DATA_DIR$INPUT_PREFIX --split-par 'hg19' --geno $Geno_Thresh --chr $CHR --export vcf --out $NO_MISSING_OUT".chr"$CHR
	#Run Phasing
	$SOFTWARE_DIR"shapeit" -V $NO_MISSING_OUT".chr"$CHR".vcf" \
		-M $REF_DIR"genetic_map_chr"$CHR"_combined_b37.txt" \
		-O $NO_MISSING_OUT".chr"$CHR".phased" \
		--thread $N_THREAD \
		--force
	$SOFTWARE_DIR"shapeit" -convert \
		--input-haps $NO_MISSING_OUT".chr"$CHR".phased" \
		--output-vcf $PHASED_DIR$INPUT_PREFIX".nomissing.chr"$CHR".phased.vcf"

	bgzip $PHASED_DIR$INPUT_PREFIX".nomissing.chr"$CHR".phased.vcf"
	#Create m3vcfs files
	$SOFTWARE_DIR"Minimac3-omp" --cpus 1 --refHap $PHASED_DIR$INPUT_PREFIX".nomissing.chr"$CHR".phased.vcf.gz" \
		--processReference --prefix $OUTPUT_DIR"m3vcfs/chr"$CHR
fi
