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
#qsub 4_split_QC2.job -v myinput=/gpfs/home/raqueld/mapping_MESA/mesa_genotypes-black.lifted_NCBI36_to_GRCh37.GH.bed,myoutdir=/stsi/raqueld/N_tests,hwe='',geno=0.1,mind=0.1 -N 4_N_mesa_genotypes-black
#qsub 4_split_QC2.job -v myinput=/stsi/raqueld/3_ancestry/6800_JHS_all_chr_sampleID_c1/6800_JHS_all_chr_sampleID_c1.lifted_hg19_to_GRCh37.GH.ancestry-5.bed,myoutdir=/stsi/raqueld/4_split_QC2,hwe='',geno=0.1,mind=0.1 -N 4_6800_JHS_all_chr_sampleID_c1

module load samtools
# module load plink2 # moved plink to required_tools

#Missingness per individual --mind N
#Missingness per marker --geno N
#Hardy-Weinberg equilibrium --hwe N

starttime=$(date +%s)

export inprefix=$(basename $myinput | sed -e 's/\.bed$//g')
export indir=$(dirname $myinput)

export plink="$HOME/required_tools/plink"
export plink2="$HOME/required_tools/plink2"

outsubdir=$(basename $myinput | sed -e 's~\.lifted.*~~g')

if [ ! -d $myoutdir/$outsubdir ]; then
    mkdir -p $myoutdir/$outsubdir
fi

cd $myoutdir/$outsubdir

if [ -z "$hwe" ]; then 
    echo "No HWE filtering applied.";
    export hweflag='' 
else 
    export hweflag=$(echo "--hwe $hwe" | tr -d '\n'); 
fi
if [ -z "$mind" ]; then 
    echo "No missingness per individual filtering applied.";
    export mindflag='' 
else 
    export mindflag=$(echo "--mind $mind" | tr -d '\n'); 
fi
if [ -z "$geno" ]; then 
    echo "No missingness per marker filtering applied.";
    export genoflag='' 
else 
    export genoflag=$(echo "--geno $geno" | tr -d '\n'); 
fi


# Shaun: retire the old split/plink-filtering cuz the inputs were split already
# #export SHELL=$(type -p bash)
# PLINK_FUN() {

#     if [ "$1" -eq 23 ]; then
#         # $plink --bfile $indir/$inprefix --chr $1,X --make-bed $hweflag $mindflag $genoflag --a1-allele  $indir/$inprefix.bim 5 2 --out $inprefix.chr$1
#         $plink2 --bfile $indir/$inprefix --chr $1,X --make-bed $hweflag $mindflag $genoflag --out $inprefix.chr$1
#     else
#         # $plink --bfile $indir/$inprefix --chr $1 --make-bed $hweflag $mindflag $genoflag --a1-allele  $indir/$inprefix.bim 5 2 --out $inprefix.chr$1
#         $plink2 --bfile $indir/$inprefix --chr $1 --make-bed $hweflag $mindflag $genoflag --out $inprefix.chr$1
#     fi

# }
# export -f PLINK_FUN


PLINK2_FUN() {
    $plink2 --bfile $indir/$inprefix.chr$1 --make-bed $hweflag $mindflag $genoflag --out $inprefix.chr$1.tmp
}
export -f PLINK2_FUN

echo "Missingness/HWE-filtering with $hweflag $mindflag $genoflag"
plinkstarttime=$(date +%s)
    parallel -j 16 PLINK2_FUN ::: {1..22}
plinkendtime=$(date +%s)


REMOVE_FUN() {
    $plink2 --bfile $inprefix.chr$1.tmp --remove $inprefix.tmp.mindrem.id --make-bed --out $inprefix.chr$1
}
export -f REMOVE_FUN

echo "Drop the samples with incomplete autosome, failed --mind (missingness of individual) for certain chromosome"
cat $inprefix.chr*.tmp.mindrem.id > $inprefix.tmp.mindrem.id
cat $inprefix.tmp.mindrem.id | sort -r | uniq > $inprefix.tmp.unique.mindrem.id
parallel REMOVE_FUN ::: {1..22}

# plink2 doesn't support merge yet, though we didn't need it necessarily here.
# echo "Merging all chromosomes..."
# #FIXING STUPID BUG WITH keep-allele-order
# for i in {2..23}; do
#     if [ -f $inprefix.chr$i.bed ]; then
#         echo $inprefix.chr$i >> $inprefix.chrlist
#     fi 
# done 
# $plink --bfile $inprefix.chr1 --merge-list $inprefix.chrlist --remove $inprefix.irem --allow-extra-chr --a1-allele $indir/$inprefix.bim 5 2 --biallelic-only --set-missing-var-ids @:#\$1:\$2 --make-bed --out $inprefix.m

if [ -f $inprefix.chrlist ]; then 
   rm $inprefix.chrlist
fi


# Cuz the next step was actually not using the merged file...
echo "Ready for next step, results ready for phasing are in $myoutdir/$outsubdir"

endtime=$(date +%s)

plinkruntime=$((plinkendtime-plinkstarttime))
overallruntime=$((endtime-starttime))

echo "plink run time: $plinkruntime"
echo "Overall run time: $overallruntime"
