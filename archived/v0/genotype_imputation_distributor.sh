#!/usr/bin/env bash

# Parameters copied to new array, preserving elements with spaces.
# Not required but defensive. BasNote "${@}" automatically drops $0.
args=( "${@}" )

# Prevent calling this function twice in the same sub shell.
if [ ! -z ${_get_opts_called_} ]
then
    echo "Can't call get_opts twice!"
    return 1
fi
_get_opts_called_=true

# Declare the array to build the translated arguments into.
declare -a rawOptions
pos=0
for elem in "${args[@]}"     # Array as list
do
    case "${elem}" in
        ### START Long options:
        "--confirm")    rawOptions[$pos]="-c" ;;
        "--wgs")    rawOptions[$pos]="-w" ;;
        "--start")   rawOptions[$pos]="-s" ;;
        "--end")   rawOptions[$pos]="-n" ;;
        "--vcf")   rawOptions[$pos]="-v" ;;
        "--out")   rawOptions[$pos]="-o" ;;
        "--ref")   rawOptions[$pos]="-r" ;;
        "--temp")    rawOptions[$pos]="-t" ;;
        ### END Long options:

        "--")   rawOptions=( "${rawOptions[@]}" "${args[@]:${pos}}" )
                break
                ;;
        "--"*)  echo "Unknown long option ${elem}"
                exit
                ;;
        *)      rawOptions[$pos]="${elem}" ;;
    esac
    pos=$(( $pos + 1 ))
done

# DEBUG can be set before calling to see the processed list of long
# options, if something is not working.
if [ ! -z ${DEBUG} ]
then
    echo "Option list after long to short processing:"
    echo "${rawOptions[@]}"
fi

while getopts "cws:n:v:o:r:t::" option "${rawOptions[@]}"
### END Option Descriptor
do
    case "$option" in

        ### START: Option variables
        "c")  confirm=true  ;;
        "w")  wgs=true  ;;
        "s")  start_from=${OPTARG}  ;;
        "n")  stop_after=${OPTARG}  ;;
        "v")  myinput="${OPTARG}" ;;
        "o")  outroot="${OPTARG}" ;;
        "r")  ref="${OPTARG}" ;;
        "t")  temp="${OPTARG}" ;;
        ### END: Option variables

        \?) echo "Unknown option -${OPTARG}"
            exit
            ;;
        \:) echo "Option -${OPTARG} requires a value."
            exit
            ;;
    esac
done

### START: Arguments variable
opt_args=( "${opt_args[@]}" "${rawOptions[@]:(( ${OPTIND} - 1 ))}" )
### STOP: Arguments variable

# Call to get_opts, with -- defining the end of options, allowing for option-like parameters




echo
echo -e " #################################### "
echo -e " ##                                ## "
echo -e " ##    Imputation / QC Pipeline    ## "
echo -e " ##          Torkamani Lab         ## "
echo -e " ##                                ## "
echo -e " ##         Author: Raquel Dias    ## "
echo -e " ##                 Shaun Chen     ## "
echo -e " ##  Last modified: 12/27/19       ## "
echo -e " ##                                ## "
echo -e " #################################### "
echo 
echo -e "Usage:    bash script.sh --vcf --out --ref --start --end (--confirm) > LOG"
echo
echo -e "          script.sh      This script"
echo -e "          --vcf -v [STR]      Full path of the input vcf file to be QCed/imputed"
echo -e "          --out -o [STR]      Path of the directory where all the output folders will be created"
echo -e "          --ref -r [STR]      Imputation reference panel (HRC or 1000G)"
echo -e "          --start -s [INT]    First step"
echo -e "          --end -e [INT]      Last step"
echo -e "          --temp -t [STR]     (Optional) Enable larger temp storage in step2 (or use PBSTMPDIR scratch folder)"
echo -e "          --wgs -w            (Optional) Enable variant down-sampling in step3 for WGS/imputed data"
echo -e "          --confirm -c        (Optional) Initiate working mode"
echo -e "          LOG                 (Optional) Log report file name"
echo
echo -e "Debug Example:  bash script.sh --vcf /mnt/stsi/stsi0/raqueld/vcf/SHARE_MESA_c2_flipfix.vcf --out /mnt/stsi/stsi0/raqueld --ref HRC --start 0 --end 1 > MESA_jobs_c1_0-1.txt"
echo -e "Working Example:  bash script.sh --vcf /mnt/stsi/stsi0/raqueld/vcf/c1_ARIC_WGS_Freeze3.vcf --out /mnt/stsi/stsi0/sfchen/dbgap/aric_wgs --ref HRC --start 2 --end 6 --temp /mnt/stsi/stsi0/sfchen/temp --wgs --confirm > ARIC_WGS_2-6.txt"
echo 





# Default settings
if [ -z $myinput ]; then echo "WARNING:  Please provide the path of input vcf file."; echo ; exit; fi
if [ -z $outroot ]; then outroot="$PWD/IMP_QC"; echo "INFO:     No OUT_ROOT, assign default path: " $outroot; fi
if [ "$ref" != "HRC" ] && [ "$ref" != "1KG" ]; then echo "WARNING:  Invalid reference panel, please use current supported HRC or 1KG."; echo; exit; fi
if [ -z $start_from ]; then start_from=999; fi
if [ -z $stop_after ]; then stop_after=999; fi
if [ -z $temp ]; then echo "INFO:     No CUSTOM_TEMP, use system scratch folder (limited to 30GB temporary storage)"; fi
if [ "$wgs" == true ]; then wgs_mode='yes'; else wgs_mode='no'; fi
if [ "$confirm" == true ]; then run=1; else run=0; fi


echo 
echo "---------------------"
echo "## Check arguments ##"
echo "---------------------"
echo "User input:  VCF_PATH: ${myinput}"
echo "             OUT_ROOT: ${outroot}"
echo "             CUSTOM_TEMP: ${temp}"
echo "User option: REF=${ref}, START=${start_from}, STOP=${stop_after}, WGS=${wgs_mode}, RUN=${run}"
echo 


###################
## Set Varaibles ##
###################
# Parsing system argument
indir=$(dirname ${myinput})
infile=$(basename ${myinput})
prefix=$(echo $infile | sed -e 's/\.vcf$//g')

# Assign submission command line
job0='qsub 0_check_vcf_build.job -v myinput=${myinput},myoutput=${outroot}/0_check_vcf_build/${prefix}.BuildChecked,copyoutput=no,gz=no -N 0_${prefix}'
job1='qsub 1_lift_vcfs_to_GRCh37.job -v myinput=${myinput},buildcheck=${outroot}/0_check_vcf_build/${prefix}.BuildChecked,myoutdir=${outroot}/1_lift,custom_temp=${temp},copyoutput=no -N 1_${prefix}'
job2='qsub 2_Genotype_Harmonizer_QC1.job -v myinput=${outroot}/1_lift/${prefix}.${lifted_code}.bed,myoutdir=${outroot}/2_GH -N 2_${prefix}'
job3='qsub 3_ancestry_analysis.job -v myinput=${outroot}/2_GH/${prefix}.${lifted_code}.GH.vcf.gz,myoutdir=${outroot}/3_ancestry,WGS=${wgs_mode} -N 3_${prefix}'
job4='qsub 4_split_QC2.job -v myinput=${outroot}/3_ancestry/${prefix}/${prefix}.${lifted_code}.GH.ancestry-${anc}.bed,myoutdir=${outroot}/4_split_QC2,geno=0.1,mind=0.05 -N 4_${prefix}'
job5='qsub 5_phase.job -v myinput=${outroot}/4_split_QC2/${prefix}/${prefix}.${lifted_code}.GH.ancestry-${anc}.chr${chrom}.bed,myoutdir=${outroot}/5_phase,reftype=${ref} -N 5_${prefix}'
job6='qsub 6_impute.job -v myinput=${outroot}/5_phase/${prefix}/${prefix}.${lifted_code}.GH.ancestry-${anc}.chr${chrom}.phased.vcf.gz,myoutdir=${outroot}/6_impute_${ref},reftype=${ref} -N 6_${prefix}'

echo "--------------------------"
echo "## Preview command line ##"
echo "--------------------------"
echo "~$ "$job0
echo "~$ "$job1
echo "~$ "$job2
echo "~$ "$job3
echo "~$ "$job4
echo "~$ "$job5
echo "~$ "$job6
echo

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
            echo ${job_ID}
                        
            sleep 0.5
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

echo "--------------------"
echo "## Job submission ##"
echo "--------------------"

job "$job0" "$run" 0 "$start_from" "$stop_after"; echo

job "$job1" "$run" 1 "$start_from" "$stop_after"; echo



{ # try

    # Get lifted status
    lifted_code=$(ls "${outroot}/1_lift/" | grep ${prefix} | grep 'lifted' | head -1 | tr '.' '\n' | grep 'lifted') &&

    {
    if [ -z $lifted_code ]; then 
        echo
        echo "WARNING: Lifted status is not ready yet, please re-submit step2-6 after step0-1 completed."
        echo
        exit
    else
        echo
        echo "INFO: Lifted status: $lifted_code"; echo

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
    fi
    }

} || { # catch
    # echo for exception 
    echo
    echo "WARNING: Lifted status is not ready yet, please re-submit step2-6 after step0-1 completed."
    echo
    exit
}

