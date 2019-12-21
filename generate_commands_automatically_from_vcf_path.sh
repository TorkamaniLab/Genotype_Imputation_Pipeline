
#example run: bash generate_commands_from_vcf_path.sh /mnt/stsi/stsi0/raqueld/vcf/SHARE_MESA_c1_flipfix.vcf /mnt/stsi/stsi0/raqueld > MESA_jobs_c2.txt
#example run: bash generate_commands_from_vcf_path.sh /mnt/stsi/stsi0/raqueld/vcf/SHARE_MESA_c2_flipfix.vcf > MESA_jobs_c1.txt
#example run: bash generate_commands_from_vcf_path.sh /mnt/stsi/stsi1/raqueld/vcf/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_flipfix.vcf > ARIC_jobs.txt


if [ -z $1 ]; then
    echo -e "Usage:    bash script.sh VCF_PATH OUT_ROOT > LOG"
    echo
    echo -e "          script.sh    This script"
    echo -e "          VCF_PATH     Full path of the input vcf file to be QCed/imputed"
    echo -e "          OUT_ROOT     Path of the directory where all the output folders will be created"
    echo -e "          LOG          Log report file name"
    echo
    echo -e "Example:  bash generate_commands_from_vcf_path.sh /mnt/stsi/stsi0/raqueld/vcf/SHARE_MESA_c2_flipfix.vcf /mnt/stsi/stsi0/raqueld > MESA_jobs_c1.txt"
    echo 
    exit
fi

if [ -z $2 ]; then
    outroot="/mnt/stsi/stsi0/raqueld"
else
    outroot=$2
fi

indir=$(dirname $1)
infile=$(basename $1)
prefix=$(echo $infile | sed -e 's/\.vcf$//g')



run=0 #whether really submit the jobs or just print out the job commands so you can debug first (recommended to debug first with run=0, then when everything is working, then run with run=1)
start_from=0 #first step to start running
stop_after=4 #last step before stop running

counter=0

echo "qsub 0_check_vcf_build.job -v myinput=$1,myoutput=$outroot/0_check_vcf_build/$prefix.BuildChecked,copyoutput=no,gz=no -N 0_$prefix"
if [ $run -eq 1 ] && [ $counter -ge $start_from ]; then
    job0=$(qsub 0_check_vcf_build.job -v myinput=$1,myoutput=$outroot/0_check_vcf_build/$prefix.BuildChecked,copyoutput=no,gz=no -N 0_$prefix)
    echo $job0
    sleep 1
fi

if [ $stop_after -eq $counter ]; then
    exit
fi

counter=$((counter+1))

if [ $counter -gt $start_from ]; then
   flag="-W depend=afterany:$job0"
else
   flag=""
fi

echo "qsub 1_lift_vcfs_to_GRCh37.job -v myinput=$1,buildcheck=$outroot/0_check_vcf_build/$prefix.BuildChecked,myoutdir=$outroot/1_lift,copyoutput=no -N 1_$prefix"
if [ $run -eq 1 ] && [ $counter -ge $start_from ]; then
    job1=$(qsub $flag 1_lift_vcfs_to_GRCh37.job -v myinput=$1,buildcheck=$outroot/0_check_vcf_build/$prefix.BuildChecked,myoutdir=$outroot/1_lift,copyoutput=no -N 1_$prefix)
    echo $job1
    sleep 1
fi

if [ $stop_after -eq $counter ]; then
    exit
fi

counter=$((counter+1))

if [ $counter -gt $start_from ]; then
   flag="-W depend=afterany:$job1"
else
   flag=""
fi

echo "qsub 2_Genotype_Harmonizer_QC1.job -v myinput=$outroot/1_lift/$prefix.lifted_NCBI36_to_GRCh37.bed,myoutdir=$outroot/2_GH -N 2_$prefix"
if [ $run -eq 1 ] && [ $counter -ge $start_from ]; then
    job2=$(qsub $flag 2_Genotype_Harmonizer_QC1.job -v myinput=$outroot/1_lift/$prefix.lifted_NCBI36_to_GRCh37.bed,myoutdir=$outroot/2_GH -N 2_$prefix)
    echo $job2
    sleep 1
fi

if [ $stop_after -eq $counter ]; then
    exit
fi

counter=$((counter+1))

if [ $counter -gt $start_from ]; then
   flag="-W depend=afterany:$job2"
else
   flag=""
fi

echo "qsub 3_ancestry_analysis.job -v myinput=$outroot/2_GH/$prefix.lifted_NCBI36_to_GRCh37.GH.vcf.gz,myoutdir=$outroot/3_ancestry -N 3_$prefix"
if [ $run -eq 1 ] && [ $counter -ge $start_from ]; then
    job3=$(qsub $flag 3_ancestry_analysis.job -v myinput=$outroot/2_GH/$prefix.lifted_NCBI36_to_GRCh37.GH.vcf.gz,myoutdir=$outroot/3_ancestry -N 3_$prefix)
    echo $job3
    sleep 1
fi

if [ $stop_after -eq $counter ]; then
    exit
fi

counter=$((counter+1))

if [ $counter -gt $start_from ]; then
   flag="-W depend=afterany:$job3"
   flag4="-W depend=afterany:$job3"
else
   flag=""
   flag4=""
fi

if [ $stop_after -eq 4 ]; then

    for i in $(echo 1 2 3 4 5 mixed); do
        echo "qsub 4_split_QC2.job -v myinput=$outroot/3_ancestry/$prefix/$prefix.lifted_NCBI36_to_GRCh37.GH.ancestry-$i.bed,myoutdir=$outroot/4_split_QC2,geno=0.1,mind=0.05 -N 4_$prefix"

        if [ $run -eq 1 ] && [ $counter -ge $start_from ]; then

           job4=$(qsub $flag4 4_split_QC2.job -v myinput=$outroot/3_ancestry/$prefix/$prefix.lifted_NCBI36_to_GRCh37.GH.ancestry-$i.bed,myoutdir=$outroot/4_split_QC2,geno=0.1,mind=0.05 -N 4_$prefix)
           echo $job4
           sleep 1
        fi
 
    done

    exit
fi

for i in $(echo 1 2 3 4 5 mixed); do
    echo "qsub 4_split_QC2.job -v myinput=$outroot/3_ancestry/$prefix/$prefix.lifted_NCBI36_to_GRCh37.GH.ancestry-$i.bed,myoutdir=$outroot/4_split_QC2,geno=0.1,mind=0.05 -N 4_$prefix"

    if [ $run -eq 1 ] && [ $counter -ge $start_from ]; then

       job4=$(qsub $flag4 4_split_QC2.job -v myinput=$outroot/3_ancestry/$prefix/$prefix.lifted_NCBI36_to_GRCh37.GH.ancestry-$i.bed,myoutdir=$outroot/4_split_QC2,geno=0.1,mind=0.05 -N 4_$prefix)
       echo $job4
       sleep 1

       counter=$((counter+1))

       if [ $counter -gt $start_from ]; then
          flag="-W depend=afterany:$job4"
       else
          flag=""
       fi

       for j in $(seq 1 1 22); do

           echo "qsub 5_phase.job -v myinput=$outroot/4_split_QC2/$prefix/$prefix.lifted_NCBI36_to_GRCh37.GH.ancestry-$i.chr$j.bed,myoutdir=$outroot/5_phase,reftype=HRC -N 5_$prefix"
           echo "qsub 6_impute.job -v myinput=$outroot/5_phase/$prefix/$prefix.lifted_NCBI36_to_GRCh37.GH.ancestry-$i.chr$j.phased.vcf.gz,myoutdir=$outroot/6_impute,reftype=HRC -N 6_$prefix"

           if [ $run -eq 1 ]; then

               job5=$(qsub $flag 5_phase.job -v myinput=$outroot/4_split_QC2/$prefix/$prefix.lifted_NCBI36_to_GRCh37.GH.ancestry-$i.chr$j.bed,myoutdir=$outroot/5_phase,reftype=HRC -N 5_$prefix)
               sleep 1

               job6=$(qsub -W depend=afterany:$job5 6_impute.job -v myinput=$outroot/5_phase/$prefix/$prefix.lifted_NCBI36_to_GRCh37.GH.ancestry-$i.chr$j.phased.vcf.gz,myoutdir=$outroot/6_impute/${prefix}/ancestry_${i},reftype=HRC -N 6_$prefix)
               sleep 1

           fi

       done
    fi
done

