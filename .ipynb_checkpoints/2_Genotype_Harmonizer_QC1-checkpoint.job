#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l mem=120gb
#PBS -q stsi 
#PBS -l walltime=540:00:00
#PBS -j oe

date

echo "Running on node:"
hostname
pwd

#myoutdir example:
#/stsi/raqueld/2_GH
#myinput example:
#qsub 2_Genotype_Harmonizer.job -v myinput=/gpfs/home/raqueld/mapping_MESA/mesa_genotypes-black.lifted_NCBI36_to_GRCh37.bed,myoutdir=/gpfs/home/raqueld/mapping_MESA,ref_path=/mnt/stsi/stsi0/raqueld/1000G -N 2_N_GH.mesa_genotypes-black

#Step 2
#Split genotype file by chromosome
#Filter SNPs by 90% call rate (per variant, not subject)
#Fix ambiguous genotype swaps
#Fix non ambiguous genotype swaps
#Remove duplicate SNPs

module load vcftools
module load samtools/1.9
module load plink2
#TODO $ which plink2

export GH=$PBS_O_WORKDIR/required_tools/GenotypeHarmonizer/GenotypeHarmonizer.jar

starttime=$(date +%s)

export inprefix=$(basename $myinput | sed -e 's/\.bed$//g')
export indir=$(dirname $myinput)

#change this path to your own reference path
if [ -z $ref_path ]; then
    export ref_path=/mnt/stsi/stsi0/raqueld/1000G
    echo "No reference path provided, check step 2 instructions in README. Using the following path as default: $ref_path"
fi

if [ ! -d $myoutdir ]; then
        mkdir -p $myoutdir
fi

cd $myoutdir

export outname=$(basename $myinput).GH


############################
# Shaun: shift split chromosome to step 1.
# #export SHELL=$(type -p bash)
# PLINK_FUN() {
#     if [ "$1" -eq 23 ]; then
#         # plink --bfile $indir/$inprefix --biallelic-only --a1-allele $indir/$inprefix.bim 5 2 --set-missing-var-ids @:#\$1:\$2 --allow-extra-chr --chr $1,X --make-bed --out $inprefix.chr$1
#         plink2 --bfile $indir/$inprefix --max-alleles 2 --set-missing-var-ids @:#\$1:\$2 --allow-extra-chr --chr $1,X --make-bed --out $inprefix.chr$1
#     else
#         # plink --bfile $indir/$inprefix --biallelic-only --a1-allele $indir/$inprefix.bim 5 2 --set-missing-var-ids @:#\$1:\$2 --allow-extra-chr --chr $1 --make-bed --out $inprefix.chr$1
#         plink2 --bfile $indir/$inprefix --max-alleles 2 --set-missing-var-ids @:#\$1:\$2 --allow-extra-chr --chr $1 --make-bed --out $inprefix.chr$1
#     fi
# }
# export -f PLINK_FUN
############################


GH_FUN () {
    if [ "$1" -eq 23 ]; then
        java -Xmx16g -jar $GH --keep --input $myinput.chr$1 \
        --ref $ref_path/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz \
        --inputType PLINK_BED --callRateFilter 0.90 --output ./$outname.chr$1

    else
        java -Xmx16g -jar $GH --keep --input $myinput.chr$1 \
        --ref  $ref_path/ALL.chr$1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
        --inputType PLINK_BED --callRateFilter 0.90 --output ./$outname.chr$1
    fi
}

export -f GH_FUN


############################
# Shaun: shift split chromosome to step 1.
# echo "Extracting chromosomes..."
# plinkstarttime=$(date +%s)
#     parallel -j 8 PLINK_FUN {1} ::: {1..23}
# plinkendtime=$(date +%s)
# echo "Done extracting chromosomes."
############################


echo "Running GH..."
GHstarttime=$(date +%s)
    parallel -j 8 GH_FUN {1} ::: {1..23}
GHendtime=$(date +%s)
echo "GH done."


############################
# Shaun: upgrade to parallel
# if [ -f $inprefix.chrlist ]; then 
#     rm $inprefix.chrlist
# fi

# mergestarttime=$(date +%s)

# for i in {2..23}; do

#     if [ -f $outname.chr$i.bed ]; then
#         echo $outname.chr$i >> $inprefix.chrlist
#         #echo $outname.chr$i.bed >> $inprefix.chrlist
#         #echo $outname.chr$i.bim >> $inprefix.chrlist
#         #echo $outname.chr$i.fam >> $inprefix.chrlist
#     fi

# done

# echo "Merging all chromosomes..."
# #FIXING STUPID BUG WITH  --keep-allele-order
# plink --bfile $outname.chr1 --merge-list $inprefix.chrlist --allow-extra-chr --a1-allele $indir/$inprefix.bim 5 2  --biallelic-only --set-missing-var-ids @:#\$1:\$2 --make-bed --out ./$outname.m
# # TODO: update to plink2 when it supported merging function.   

# #RR moving fix ref and introducing remove duplicates
# #DISCOVERED THAT --keep-allele-order IS NOT WORKING, DISABLING ALLELE ORDER CHANGE BY FORCE!!!!!!
# # plink --bfile $outname.m --recode vcf bgz --a1-allele $outname.m.bim 5 2 --double-id --set-missing-var-ids @:#\$1:\$2 --out $outname.0
# plink2 --bfile $outname.m --export vcf-4.2 bgz --double-id --set-missing-var-ids @:#\$1:\$2 --out $outname.0

# tabix -p vcf $outname.0.vcf.gz
# #rm $outname.m.bed $outname.m.bim $outname.m.fam

# #echo "Doing last check for allele swaps"
# bcftools +fixref $outname.0.vcf.gz --threads 16 -Oz -o $outname.1.vcf.gz -- -f $ref_path/human_g1k_v37.fasta -m flip -d

# tabix -p vcf $outname.1.vcf.gz
# #rm $outname.0.vcf.gz

# bcftools view -Ou -c 2 $outname.1.vcf.gz | bcftools norm -m -any | bcftools norm --threads 16 -Oz -o $outname.vcf.gz -d both -f $ref_path/human_g1k_v37.fasta

# tabix -p vcf $outname.vcf.gz

# #rm $outname.1.vcf.gz
# #rm $outname.0.vcf.gz
############################


FIXREF_FUN () {
    echo "Merging all chromosomes..."
    # NOTE: --id-delim can no longer be used with --const-fid or --double-id.
    plink2 --bfile $outname.chr$1 --max-alleles 2 --set-missing-var-ids @:#\$1:\$2 --export vcf-4.2 bgz --double-id --out $outname.chr$1.0

    tabix -p vcf $outname.chr$1.0.vcf.gz
    #rm $outname.m.bed $outname.m.bim $outname.m.fam

    echo "Doing last check for allele swaps"
    bcftools +fixref $outname.chr$1.0.vcf.gz --threads 16 -Oz -o $outname.chr$1.1.vcf.gz -- -f $ref_path/human_g1k_v37.fasta -m flip -d

    tabix -p vcf $outname.chr$1.1.vcf.gz
    #rm $outname.0.vcf.gz

    bcftools view -Ou -c 2 $outname.chr$1.1.vcf.gz | bcftools norm -m -any | bcftools norm --threads 16 -Oz -o $outname.chr$1.vcf.gz -d both -f $ref_path/human_g1k_v37.fasta

    tabix -p vcf $outname.chr$1.vcf.gz

    #rm $outname.1.vcf.gz
    #rm $outname.0.vcf.gz
}
export -f FIXREF_FUN

echo "Running FIXREF..."
fixrefstarttime=$(date +%s)
    parallel FIXREF_FUN {1} ::: {1..22}
fixrefendtime=$(date +%s)
echo "FIXREF done."


CLEAN_FUN () {
if [ -f $outname.chr$1.bed ]; then
    rm $outname.chr$1.bed
    rm $outname.chr$1.bim
    rm $outname.chr$1.fam
fi
if [ -f $inprefix.chr$1.bed ]; then
    rm $inprefix.chr$1.*
fi
}
export -f CLEAN_FUN

echo "Cleaning temporary files..."
    parallel CLEAN_FUN {1} ::: {1..22}
echo "Cleaning done."


endtime=$(date +%s)

# plinkruntime=$((plinkendtime-plinkstarttime))
GHruntime=$((GHendtime-GHstarttime))
# mergeruntime=$((mergeendtime-mergestarttime))
fixreftime=$((fixrefendtime-fixrefstarttime))
overallruntime=$((endtime-starttime))

# echo "plink run time: $plinkruntime"
echo "GH run time: $GHruntime"
# echo "merge run time: $mergeruntime"
echo "fixref run time: $fixreftime"
echo "Overall run time: $overallruntime"

#output example
#ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out.GH.chr1.vcf.gz