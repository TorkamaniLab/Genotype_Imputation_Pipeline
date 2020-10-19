#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=64gb
# #PBS -q stsi
#PBS -l walltime=340:00:00
#PBS -j oe

echo "Running on node:"
hostname

module load R/3.5.1

#myinput example: /stsi/raqueld/vcf/6800_JHS_all_chr_sampleID_c2.vcf
#myoutput example: /stsi/raqueld/0_check_vcf_build/6800_JHS_all_chr_sampleID_c2.BuildChecked
#How to run Example
#qsub 0_check_vcf_build.job -v myinput=/stsi/raqueld/vcf/6800_JHS_all_chr_sampleID_c2.vcf,myoutput=/stsi/raqueld/0_check_vcf_build/6800_JHS_all_chr_sampleID_c2.BuildChecked,gz=yes -N 0_6800_JHS_all_chr_sampleID_c2

# READ_ARR() {
#     # create filepath array when input txt file list.
#     declare -A my_arr
#     while read line; do 
#         chrom=$(echo $line | awk '{print$1}')
#         fp=$(echo $line | awk '{print$2}')
#     #     echo $chrom and $fp
#         my_arr[$chrom]=${fp}
#     done < $1
# }
# export -f READ_ARR


# Get only filename of chr1 to check the genome build if input vcf already split
if [[ ${myinput} == *.txt ]]; then
    echo "Input autosome list detected, run parallel pipeline."
    my_chr1=$(cat ${myinput} | awk '$1==1 {print$2}')
    mydirname=$(dirname ${my_chr1})
    filename=$(basename ${my_chr1})
    echo $mydirname/$filename
elif [[ ${myinput} == *.vcf.gz ]]; then
    echo "Merged vcf.gz file detected, chr1 will be extract in the Rscript."
    mydirname=$(dirname $myinput)
    filename=$(basename $myinput)    
else
    echo "Invalid file format, please check the input."
    exit
fi

echo "Running vcf build check"

# required_tools folder need to locate with this job script
check_vcf_build=$PBS_O_WORKDIR/required_tools/check_vcf_build/check_vcf_build.R
if [ ! -f $check_vcf_build ]; then
    echo "check_vcf_build.R not found at ${check_vcf_build}. Make sure you downloaded check_vcf_build.R and changed the check_vcf_build variable inside this job script, so it matches to the path where your check_vcf_build.R if located."
    exit
fi

#check if output directory exists
myoutdir=$(dirname $myoutput)
if [ ! -d $myoutdir ]; then
    mkdir -p $myoutdir
fi

# Shaun: disable decompressed input here anymore.
# 1) Enable `less` in extracting positions in the R script.
# gzstarttime=$(date +%s)
# if [ "$gz" == "yes" ]; then
#     echo "Extracting $filename to ${filename/.gz/}"
#     cd $mydirname
#         bgzip -dc -@ 16 $filename > ${filename/.gz/}
#     filename=${filename/.gz/}
# fi
# gzendtime=$(date +%s)

vcfstarttime=$(date +%s)
Rscript $check_vcf_build $mydirname/$filename > $myoutput
vcfendtime=$(date +%s)

vcfruntime=$((vcfendtime-vcfstarttime))
# gzruntime=$((gzendtime-gzstarttime))

echo "vcf check run time: $vcfruntime"
# echo "Extract run time: $gzruntime"

