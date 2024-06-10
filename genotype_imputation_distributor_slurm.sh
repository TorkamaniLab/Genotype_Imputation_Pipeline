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
        "--hwe")    rawOptions[$pos]="-h" ;;       
        "--archived")    rawOptions[$pos]="-a" ;;
        "--start")   rawOptions[$pos]="-s" ;;
        "--end")   rawOptions[$pos]="-n" ;;
        "--vcf")   rawOptions[$pos]="-v" ;;
        "--out")   rawOptions[$pos]="-o" ;;
        "--ref")   rawOptions[$pos]="-r" ;;
        "--lai-ref")   rawOptions[$pos]="-l" ;;
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

while getopts "cwha:s:n:v:o:r:l:t::" option "${rawOptions[@]}"
### END Option Descriptor
do
    case "$option" in

        ### START: Option variables
        "c")  confirm_run=true  ;;
        "w")  wgs_mode=true  ;;
        "h")  hwe_mode=true  ;;
        "a")  archived=${OPTARG}  ;;
        "s")  start_from=${OPTARG}  ;;
        "n")  stop_after=${OPTARG}  ;;
        "v")  myinput="${OPTARG}" ;;
        "o")  outroot="${OPTARG}" ;;
        "r")  ref_mode="${OPTARG}" ;;
        "l")  LAI_ref="${OPTARG}" ;;
        "t")  tempdir="${OPTARG}" ;;
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
echo -e " ##  Last modified: 01/12/21       ## "
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
echo -e "          --temp -t [STR]     (Optional) Enable larger temp storage in step2 (or use TMP scratch folder)"
echo -e "          --wgs -w            (Optional) Enable variant down-sampling in step3 for WGS/imputed data"
echo -e "          --hwe -h            (Optional) Enable Hardy–Weinberg equilibrium (HWE) filering in step4 for non-mixed population 'plink --hwe 1e-10'"
echo -e "          --confirm -c        (Optional) Initiate working mode"
# echo -e "          --lai-ref -l [STR]  (todo) Local ancestry inference reference panel (1000G, HRC-1000G, HGDP or SGDP)" #TODO
# echo -e "          --archived -a [INT] (todo) Enable archived version (start from 0); use the latest one if no flag provided"  #TODO
echo -e "          LOG                 (Optional) Log report file name"
echo
echo -e 'Prepare input:  for chrom in {1..22}; do printf "${chrom}\\t$(ls $(pwd)/[##inprefix##]*chr${chrom}.vcf.gz)\\n"; done > [##inprefix##].txt'
echo
echo -e "Debug Example:  bash script.sh --vcf /mnt/stsi/stsi0/raqueld/vcf/SHARE_MESA_c2_flipfix.vcf --out /mnt/stsi/stsi0/raqueld --ref HRC --start 0 --end 1 > MESA_jobs_c1_0-1.txt"
echo -e "Working Example:  bash script.sh --vcf /mnt/stsi/stsi0/raqueld/vcf/c1_ARIC_WGS_Freeze3.vcf --out /mnt/stsi/stsi0/sfchen/dbgap/aric_wgs --ref HRC --start 2 --end 6 --temp /mnt/stsi/stsi0/sfchen/temp --wgs --hwe --confirm > ARIC_WGS_2-6.txt"
echo 




echo 
echo "---------------------"
echo "## Check arguments ##"
echo "---------------------"
# Default settings
if [ -z $myinput ]; then echo "WARNING: Please provide the path of input vcf file."; echo ; exit; fi
if [ -z $outroot ]; then outroot="$PWD/IMP_QC"; echo "INFO: No OUT_ROOT, assign default path: " $outroot; fi
if [ "$ref_mode" != "HRC" ] && [ "$ref_mode" != "1KG" ]; then echo "WARNING: Invalid reference panel, please use current supported HRC or 1KG."; echo; exit; fi
# if [ "$LAI_ref" != "HRC" ] && [ "$LAIref" != "1KG" ]; then echo "WARNING:  Invalid reference panel, please use current supported HRC or 1KG."; echo; exit; fi
if [ -z $start_from ] || ! [[ "$start_from" =~ ^[0-9]+$ ]]; then echo "WARNING: Invalid start step number."; echo; exit; fi
if [ -z $stop_after ] || ! [[ "$stop_after" =~ ^[0-9]+$ ]]; then echo "WARNING: Invalid end step number."; echo; exit; fi
if [ $start_from -gt $stop_after ]; then echo "WARNING: end step number should be larger then start one."; echo; exit; fi
# if [ -z $archived ]; then echo "INFO: No archived: use the latest pipeline."; else archived=archived/v${archived}/; fi
if [ -z $tempdir ]; then echo "INFO: No CUSTOM_TEMP: use system scratch folder (limited to 30GB temporary storage)"; fi
if [ "$wgs_mode" == true ]; then wgs='yes'; else wgs='no'; fi
if [ "$hwe_mode" == true ]; then hwe='yes'; else hwe='no'; fi
if [ "$confirm_run" == true ]; then run=1; else run=0; fi
echo "User input:  VCF_PATH: ${myinput}"
echo "             OUT_ROOT: ${outroot}"
echo "             CUSTOM_TEMP: ${tempdir}"
echo "User option: REF=${ref_mode}, START=${start_from}, STOP=${stop_after}, WGS=${wgs}, HWE=${hwe}, RUN=${run}"
echo 


###################
## Set Varaibles ##
###################
# Parsing system argument
indir=$(dirname ${myinput})
infile=$(basename ${myinput})
prefix=$(echo $infile | sed -e 's/\.vcf.gz$//g' | sed -e 's/\.txt$//g' )

# Assign submission command line
schedular='sbatch'
job0='--export=ALL,myinput=${myinput},myoutput=${outroot}/0_check_vcf_build/${prefix}.BuildChecked --job-name=0_${prefix} --out=0_${prefix}.o%j 0_check_vcf_build.slurm.sh'
job1='--export=ALL,myinput=${myinput},buildcheck=${outroot}/0_check_vcf_build/${prefix}.BuildChecked,myoutdir=${outroot}/1_lift,custom_temp=${outroot}/1_lift/temp --job-name=1_${prefix} --out=1_${prefix}.o%j 1_lift_vcfs_to_GRCh37.slurm.sh'
job2='--export=ALL,myinput=${outroot}/1_lift/${prefix}.${lifted_code},myoutdir=${outroot}/2_GH --job-name=2_${prefix} --out=2_${prefix}.o%j 2_Genotype_Harmonizer_QC1.slurm.sh'
job3='--export=ALL,myinput=${outroot}/2_GH/${prefix}.${lifted_code}.GH,myoutdir=${outroot}/3_ancestry,WGS=${wgs} --job-name=3_${prefix} --out=3_${prefix}.o%j 3_ancestry_analysis.slurm.sh'
job4='--export=ALL,myinput=${outroot}/3_ancestry/${prefix}/${prefix}.${lifted_code}.GH.ancestry-${anc},myoutdir=${outroot}/4_split_QC2,geno=0.1,mind=0.05 --job-name=4_${prefix} --out=4_${prefix}.o%j 4_split_QC2.slurm.sh'
job4_hwe='--export=ALL,myinput=${outroot}/3_ancestry/${prefix}/${prefix}.${lifted_code}.GH.ancestry-${anc},myoutdir=${outroot}/4_split_QC2,geno=0.1,mind=0.05,hwe=1e-10 --job-name=4_${prefix} --out=4_${prefix}.o%j 4_split_QC2.slurm.sh'
job5='--export=ALL,myinput=${outroot}/4_split_QC2/${prefix}/${prefix}.${lifted_code}.GH.ancestry-${anc}.chr${chrom}.bed,myoutdir=${outroot}/5_phase,reftype=${ref_mode} --job-name=5_${prefix} --out=5_${prefix}.o%j 5_phase.slurm.sh'
job6='--export=ALL,myinput=${outroot}/5_phase/${prefix}/${prefix}.${lifted_code}.GH.ancestry-${anc}.chr${chrom}.phased.vcf.gz,myoutdir=${outroot}/6_impute_${ref_mode},reftype=${ref_mode} --job-name=6_${prefix} --out=6_${prefix}.o%j 6_impute.slurm.sh'

echo "--------------------------"
echo "## Preview command line ##"
echo "--------------------------"
echo "~\$ "$schedular $job0
echo "~\$ "$schedular $job1
echo "~\$ "$schedular $job2
echo "~\$ "$schedular $job3
echo "~\$ "$schedular $job4
echo "~\$ "$schedular $job5
echo "~\$ "$schedular $job6
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

make_job() {
    schedular="$1" job="$2" run="$3" step="$4" start_from="$5" stop_after="$6" anc="$7" chrom="$8"

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
        job_command=$(eval echo ${schedular} ${depend_flag} ${job})
        echo $job_command
        
        # run if debugging mode disabled, recording job_ID
        if [ $run -eq 1 ]; then
            job_ID=""
            
            # submit the job and save jobID
            submit_job=$(${job_command})
            echo ${submit_job}
            
            job_ID=$(echo ${submit_job} | tr ' ' '\n' | tail -1)
                        
            sleep 0.5
        fi
        
        # send jobID to flag array (empty is fine for -W)
        next_flag="--dependency=afterok:${job_ID}"
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

make_job "$schedular" "$job0" "$run" 0 "$start_from" "$stop_after"; echo

make_job "$schedular" "$job1" "$run" 1 "$start_from" "$stop_after"; echo



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

        make_job "$schedular" "$job2" "$run" 2 "$start_from" "$stop_after"; echo

        make_job "$schedular" "$job3" "$run" 3 "$start_from" "$stop_after"; echo

        if [ "$hwe" == 'yes' ]; then
            for anc in $ancestry; do
                if [ "$anc" == "mixed" ]; then
                    make_job "$schedular" "$job4" "$run" 4 "$start_from" "$stop_after" "$anc"
                else
                    make_job "$schedular" "$job4_hwe" "$run" 4 "$start_from" "$stop_after" "$anc"
                fi
            done; echo               
        else
            for anc in $ancestry; do
                make_job "$schedular" "$job4" "$run" 4 "$start_from" "$stop_after" "$anc"
            done; echo
        fi



        for anc in $ancestry; do
            for chrom in {1..22}; do
                make_job "$schedular" "$job5" "$run" 5 "$start_from" "$stop_after" "$anc" "$chrom"
            done
        done; echo

        for anc in $ancestry; do
            for chrom in {1..22}; do
                make_job "$schedular" "$job6" "$run" 6 "$start_from" "$stop_after" "$anc" "$chrom"
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

