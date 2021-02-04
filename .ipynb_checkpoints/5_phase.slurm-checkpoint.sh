#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=540:00:00
#SBATCH --mem=100G


#myoutdir example:
#/stsi/raqueld/5_phase

#job example:
#qsub 5_phase.job -v myinput=/stsi/raqueld/N_tests/aric_genotypes-black/aric_genotypes-black.lifted_NCBI36_to_GRCh37.GH.chr1.bed,myoutdir=/stsi/raqueld/5_N_tests,reftype=HRC -N 5_N_mesa_genotypes-black
#qsub 5_phase.job -v myinput=/stsi/raqueld/4_split_QC2/unique/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed.lifted_NCBI36_to_GRCh37.GH.ancestry-3.chr22.bed,myoutdir=/stsi/raqueld/5_phase,reftype=HRC -N 5_ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed
#the input must have the suffix *.lifted*.chr1.bed, *lifted*.chr2.bed, *.lifted*.chr3.bed, etc. The previous steps in the pipeline generate those suffixes automatically, but keep these suffixes in mind if you are running this step as a stand alone tools, without running the previous steps


date
echo "Running on node:"
hostname
pwd


module purge
module load samtools

export plink="$SLURM_SUBMIT_DIR/required_tools/plink"
export plink2="$SLURM_SUBMIT_DIR/required_tools/plink2"
export eagle="$SLURM_SUBMIT_DIR/required_tools/Eagle_v2.4.1/eagle"


starttime=$(date +%s)

inprefix=$(basename $myinput | sed -e 's/\.bed$//g')
indir=$(dirname $myinput)
mychr=$(echo $inprefix | sed -e 's/.*\.chr//g')
mymap="$SLURM_SUBMIT_DIR/required_tools/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz"

echo "Chromosome $mychr"
echo "Input prefix $inprefix"

outsubdir=$(basename $myinput | sed -e 's/\.lifted.*//g')

if [ ! -d $myoutdir/$outsubdir ]; then
	mkdir -p $myoutdir/$outsubdir
fi

cd $myoutdir/$outsubdir

if [ "$reftype" == "1KG" ]; then
    myref=/mnt/stsi/stsi0/raqueld/1000G/ALL.chr$mychr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
    echo "Using 1KG reference, file: $myref"
else
    myref=/mnt/stsi/stsi0/raqueld/HRC/HRC.r1-1.EGA.GRCh37.chr$mychr.haplotypes.bcf
    echo "Using HRC reference panel, file: $myref"
fi

plinkstarttime=$(date +%s)

#FIXING STUPID BROKEN KEEP-ALLELE-ORDER
# $plink --bfile  $indir/$inprefix --recode vcf bgz --a1-allele $indir/$inprefix.bim 5 2 --set-missing-var-ids @:#\$1:\$2 --out $inprefix
$plink2 --bfile  $indir/$inprefix --export vcf-4.2 bgz --set-missing-var-ids @:#\$1:\$2 --out $inprefix

tabix -f -p vcf $inprefix.vcf.gz

plinkendtime=$(date +%s)

eaglestarttime=$(date +%s)

$eagle --vcfTarget=$inprefix.vcf.gz \
--vcfRef=$myref \
--noImpMissing \
--geneticMapFile=$mymap \
--Kpbwt=100000 --numThreads=16 \
--chrom=$mychr --allowRefAltSwap \
--outPrefix=$inprefix.phased

eagleendtime=$(date +%s)

echo "Ready for next step, results ready for imputation in $myoutdir/$outsubdir"
echo "Output prefix $inprefix.phased"

endtime=$(date +%s)

eagleruntime=$((eagleendtime-eaglestarttime))
plinkruntime=$((plinkendtime-plinkstarttime))
overallruntime=$((endtime-starttime))

echo "Plink run time: $plinkruntime"
echo "Eagle run time: $eagleruntime"
echo "Overall run time: $overallruntime"
