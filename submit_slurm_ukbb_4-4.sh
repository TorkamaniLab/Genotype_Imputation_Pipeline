#!/bin/bash

parmax=$1

echo "setting max job number as ${parmax}"

for i in {022..048}; do
    
    while :; do
    
        jobcount=$(squeue -u sfchen | grep -w "shared" | wc -l)
        echo "current total job number: $jobcount - $(date +"%c")"
        
        if [ "$jobcount" -le "$parmax" ]; then
            echo "submit batch ${i}"
            bash genotype_imputation_distributor_slurm.sh \
              --vcf /mnt/stsi/stsi0/sfchen/UKBB/1_split_QC/ukb_hap_v2_${i}.txt \
              --out /mnt/stsi/stsi0/sfchen/UKBB/1_split_QC/ukbb_hap_v2_${i} \
              --ref HRC \
	      --hwe \
              --start 4 \
              --end 4 \
              --temp /mnt/stsi/stsi0/sfchen/temp --confirm > UKBB_hap_v2_${i}_4-4.txt
            break
        fi
        sleep 300

    done
done
