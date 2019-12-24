#!/usr/bin/env bash

echo
echo -e " #################################### "
echo -e " ##                                ## "
echo -e " ##    Imputation / QC Pipeline    ## "
echo -e " ##          Torkamani Lab         ## "
echo -e " ##                                ## "
echo -e " ##         Author: Raquel Dias    ## "
echo -e " ##                 Shaun Chen     ## "
echo -e " ##  Last modified: 12/24/19       ## "
echo -e " ##                                ## "
echo -e " #################################### "
echo 
echo -e "Usage:    bash script.sh VCF_PATH OUT_ROOT START STOP REF [CONFIRM] > LOG"
echo
echo -e "          script.sh      This script"
echo -e "          - VCF_PATH     Full path of the input vcf file to be QCed/imputed"
echo -e "          - OUT_ROOT     Path of the directory where all the output folders will be created"
echo -e "          - REF          Imputation reference panel (HRC or 1000G)"
echo -e "          - START        First step"
echo -e "          - STOP         Last step"
echo -e "          - [CONFIRM]    Any string to disable debugging mode"
echo -e "          - LOG          Log report file name"
echo
echo -e "Debug Example:  bash generate_commands_from_vcf_path.sh /mnt/stsi/stsi0/raqueld/vcf/SHARE_MESA_c2_flipfix.vcf /mnt/stsi/stsi0/raqueld 0 6 HRC > MESA_jobs_c1.txt"
echo 
echo -e "NOTE: When everything is working, then add CONFIRM flag to run with run = 1"
echo 

# Default settings
if [ -z $1 ]; then exit; else myinput=$1; fi
if [ -z $2 ]; then outroot="/[OUTPUT_DIRECTORY]/"; else outroot=$2; fi
if [ -z $3 ]; then ref="HRC"; else ref=$3; fi
if [ -z $4 ]; then start_from=999; else start_from=$4; fi
if [ -z $5 ]; then stop_after=999; else stop_after=$5; fi
if [ -z $6 ]; then run=0; else run=1; fi

echo "---------------"
echo "Check arguments"
echo "---------------"
echo "User input:  VCF_PATH: $1"
echo "             OUT_ROOT: $2"
echo "User option: REF=$3, START=$4, STOP=$5, RUN=${run}"
echo 


###################
## Set Varaibles ##
###################
# Parsing system argument
indir=$(dirname ${myinput})
infile=$(basename ${myinput})
prefix=$(echo $infile | sed -e 's/\.vcf$//g')

# Assign submission command line
echo "--------------------"
echo "Preview command line"
echo "--------------------"
job0='qsub 0_check_vcf_build.job -v myinput=${myinput},myoutput=${outroot}/0_check_vcf_build/${prefix}.BuildChecked,copyoutput=no,gz=no -N 0_${prefix}'
job1='qsub 1_lift_vcfs_to_GRCh37.job -v myinput=${myinput},buildcheck=${outroot}/0_check_vcf_build/${prefix}.BuildChecked,myoutdir=${outroot}/1_lift,copyoutput=no -N 1_${prefix}'
job2='qsub 2_Genotype_Harmonizer_QC1.job -v myinput=${outroot}/1_lift/${prefix}.${lifted_code}.bed,myoutdir=${outroot}/2_GH -N 2_${prefix}'
job3='qsub 3_ancestry_analysis.job -v myinput=${outroot}/2_GH/${prefix}.${lifted_code}.GH.vcf.gz,myoutdir=${outroot}/3_ancestry -N 3_${prefix}'
job4='qsub 4_split_QC2.job -v myinput=${outroot}/3_ancestry/${prefix}/${prefix}.${lifted_code}.GH.ancestry-${anc}.bed,myoutdir=${outroot}/4_split_QC2,geno=0.1,mind=0.05 -N 4_${prefix}'
job5='qsub 5_phase.job -v myinput=${outroot}/4_split_QC2/${prefix}/${prefix}.${lifted_code}.GH.ancestry-${anc}.chr${chrom}.bed,myoutdir=${outroot}/5_phase,reftype=${ref} -N 5_${prefix}'
job6='qsub 6_impute.job -v myinput=${outroot}/5_phase/${prefix}/${prefix}.${lifted_code}.GH.ancestry-${anc}.chr${chrom}.phased.vcf.gz,myoutdir=${outroot}/6_impute,reftype=${ref} -N 6_${prefix}'

# empty depend step flag by default
flag=""

# ancestry array
ancestry="1 2 3 4 5 mixed"


#####################
## Define Function ##
#####################
# declare dictionary for depend flag
declare -A flag_arr

job() {
    job="$1" run="$2" step="$3" start_from="$4" stop_after="$5" anc="$6" chrom="$7"

    # check if step isin range
    if [ $step -ge $start_from ] && [ $step -le $stop_after ]; then
        # get the dependancy job_ID from previous step
        last_step=$(expr ${step} - 1)

        # if it's step 5, depend flag need ancestry from step 4
        if [ $step -eq 5 ]; then
            depend_flag=${flag_arr["${last_step}${anc}"]}
#             echo depend_flag ${last_step}${anc}
        # if it's step 6, depend flag need ancestry/chromosome from step 5
        elif [ $step -eq 6 ]; then
            depend_flag=${flag_arr["${last_step}${anc}${chrom}"]}
#             echo depend_flag ${last_step}${anc}${chrom}
        # else, just get from previous step
        else
            depend_flag=${flag_arr["${last_step}"]}
#             echo depend_flag ${last_step}
        fi

        # recompose the qsub command line with proper variables
        job=$(eval echo ${job} ${depend_flag})
        echo $job
        
        # run if debugging mode disabled, recording job_ID
        if [ $run -eq 1 ]; then
            job_ID=""
            
            # submit the job and save jobID
            job_ID=$($job)
            echo ${job0_ID}
                        
            sleep 1
        fi
        
        # send jobID to flag array (empty is fine for -W)
        next_flag="-W depend=afterany:${job_ID}"
        flag_arr["${step}${anc}${chrom}"]=$next_flag
#         echo save_step ${step}${anc}${chrom}

    else
        if [ -z $anc ]; then anc="ALL"; fi
        if [ -z $chrom ]; then chrom="ALL"; fi
        echo "# Skip step-${step}_ancestry-${anc}_chr${chrom}"
    fi
}


#################
## Main Script ##
#################

echo "--------------"
echo "Job submission"
echo "--------------"

job "$job0" "$run" 0 "$start_from" "$stop_after"; echo

job "$job1" "$run" 1 "$start_from" "$stop_after"; echo

# Get lifted status
lifted_code=$(ls "${outroot}/1_lift/" | grep ${prefix} | grep 'lifted' | head -1 | tr '.' '\n' | grep 'lifted')
echo "# Lifted status: $lifted_code"; echo

job "$job2" "$run" 2 "$start_from" "$stop_after"; echo

job "$job3" "$run" 3 "$start_from" "$stop_after"; echo

for anc in $ancestry; do
    job "$job4" "$run" 4 "$start_from" "$stop_after" "$anc"
done; echo

for anc in $ancestry; do
    for chrom in {1..22}; do
        job "$job5" "$run" 5 "$start_from" "$stop_after" "$anc" "$chrom"
    done
done; echo

for anc in $ancestry; do
    for chrom in {1..22}; do
        job "$job6" "$run" 6 "$start_from" "$stop_after" "$anc" "$chrom"
    done
done; echo

