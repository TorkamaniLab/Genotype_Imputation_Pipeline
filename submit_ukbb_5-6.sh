#!/bin/bash

parmax=$1

echo "setting max job number as ${parmax}"

for i in {008..048}; do
    
    while :; do
    
        jobcount=$(qstat -u sfchen | grep -w "stsi" | wc -l)
        echo "current total job number: $jobcount - $(date +"%c")"
        
        if [ "$jobcount" -le "$parmax" ]; then
            echo "submit batch ${i}"
            bash genotype_imputation_distributor.sh \
              --vcf /mnt/stsi/stsi0/sfchen/UKBB/split_QC/ukb_hap_v2_${i}.txt \
              --out /mnt/stsi/stsi0/sfchen/UKBB/split_QC/ukbb_hap_v2_${i} \
              --ref HRC \
              --start 5 \
              --end 6 \
              --temp /mnt/stsi/stsi0/sfchen/temp \
              --confirm > UKBB_hap_v2_${i}_5-6.txt
            break
        fi
        sleep 300

    done
done







# #!/bin/bash
# #PBS -l nodes=1:ppn=1
# #PBS -l mem=2gb
# #PBS -q sata
# #PBS -l walltime=600:00:00
# #PBS -j oe



# #example qsub job_submitter.job -v parmax=4,joblist=eMERGE_file_list_to_reverse_imputation.txt.0000.jobs -N job.submitter0

# overallstarttime=$(date +%s)

# echo "Running on node:"
# hostname

# cd /gpfs/home/raqueld/scripts

# count=1

# while read line; do

# 	while :
# 	do

# 		jobcount=$(qstat -u raqueld | grep -w "R\|Q" | wc -l)
# 		if [ "$jobcount" -le "$parmax" ]; then
# 			echo job count is $jobcount, submitting more:
# 			echo $line
# 			sleep 3
# 			$line
# 			sleep 3
# 			echo "$count jobs submitted."
# 			count=$(expr $count + 1)
# 			break
# 		fi
# 		sleep 3;
# 	done

# done < $joblist

# overallendtime=$(date +%s)

# overallruntime=$((overallendtime-overallstarttime))

# echo "Overall run time: $overallruntime"