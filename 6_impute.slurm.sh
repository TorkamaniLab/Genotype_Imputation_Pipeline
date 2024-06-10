#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=540:00:00
#SBATCH --mem=100G


#example 
#qsub 6_impute.job -v myinput=/stsi/raqueld/5_N_tests/mesa_genotypes-white/mesa_genotypes-white.lifted_NCBI36_to_GRCh37.GH.chr18.phased.vcf.gz,myoutdir=/stsi/raqueld/6_N_tests,reftype=HRC -N 6_mesa_genotypes-white.lifted_NCBI36_to_GRCh37.GH.chr18
# if running this step as stand alone tool, input file must have the suffix .lifted_*.chr*.phased.vcf.gz, otherwise the pipeline wont work, if you use the previous step to generate this input file, then it will work fine.


date
echo "Running on node:"
hostname
pwd
echo "myinput $myinput"


module purge
module load samtools
module load minimac4/1.0.2
imputationSoftware=minimac4


inprefix=$(basename $myinput | sed -e 's/\.phased\.vcf\.gz//g')
indir=$(dirname $myinput)
mychr=$(echo $inprefix | sed -e 's/.*\.chr//g')

echo "Chromosome $mychr"
echo "Input prefix $inprefix"
subdir=$(basename $myinput | sed -e 's/\.lifted.*//g')
subanc=$(basename $myinput | tr '.' '\n' | grep "ancestry")
outsubdir=$subdir/$subanc
echo "Working directory: ${myoutdir}/${outsubdir}"

if [ ! -d $myoutdir/$outsubdir ]; then
        mkdir -p $myoutdir/$outsubdir
        # mkdir -p $myoutdir
fi

cd $myoutdir/$outsubdir
# cd $myoutdir


if [ "$reftype" == "1KG" ]; then
    myref=/mnt/stsi/stsi3/Internal/1000G/ref_panel/hg19/m3vcf_erate_rec/ALL.chr$mychr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.m3vcf.gz
    echo "Using 1KG reference panel, file: ${myref}"
else
    myref=/mnt/stsi/stsi3/Internal/HRC/ref_panel/hg19/m3vcf_erate_rec/HRC.r1-1.EGA.GRCh37.chr$mychr.haplotypes.m3vcf.gz
    echo "Using HRC reference panel, file: ${myref}"
fi

echo "Running Imputation"
$imputationSoftware --refHaps $myref \
         --haps $myinput \
         --prefix imputed_$mychr \
         --ignoreDuplicates \
         --minRatio 0.01 --ChunkLengthMb 30.00 --ChunkOverlapMb 3.00 --cpus 16

echo "Running tabix"
tabix -p vcf imputed_$mychr.dose.vcf.gz
