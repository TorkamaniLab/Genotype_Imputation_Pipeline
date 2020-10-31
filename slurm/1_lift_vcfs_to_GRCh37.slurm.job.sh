#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --partition=em
#SBATCH --time=340:00:00


# Commented out torque header
: '
#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l mem=120gb
#PBS -q stsi
#PBS -l walltime=340:00:00
#PBS -j oe
'

echo "Running on node:"
hostname

module load plink2
module load ucsc_tools/373
# module load vcftools
module load samtools # update for SLURM

#Examples of the variables needed (-v)
#myinput=/stsi/raqueld/vcf/6800_JHS_all_chr_sampleID_c2.vcf
#buildcheck=/stsi/raqueld/0_check_vcf_build/6800_JHS_all_chr_sampleID_c2.BuildChecked
#myoutdir=/stsi/raqueld/1_lift

#How to run Example
#qsub 1_lift_vcfs_to_GRCh37.job -v myinput=/stsi/raqueld/vcf/6800_JHS_all_chr_sampleID_c2.vcf,buildcheck=/stsi/raqueld/0_check_vcf_build/6800_JHS_all_chr_sampleID_c2.BuildChecked,myoutdir=/stsi/raqueld/1_lift,custom_temp=/mnt/stsi/stsi0/raqueld/tmp -N 1_6800_JHS_all_chr_sampleID_c2


export filename=$(basename $buildcheck)
export inprefix=${filename/.BuildChecked/}


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
# export -f CHECK_FORMAT


CHECK_SORT () {
    # Check if input is indexed (also check sort)
    chrom=$(echo $1 | awk '{print$1}')
    fp=$(echo $1 | awk '{print$2}')
    
    if [ ! -e ${fp}.tbi ]; then
        echo "input not tabixed"
        { # try
            echo "tabixing..."
            tabix -p vcf ${fp}
            chmod 770 ${fp}.tbi
            echo "tabix done, skip sorting."
            cp ${fp} ./${inprefix}.chr${chrom}.sorted.vcf.gz
            cp ${fp}.tbi ./${inprefix}.chr${chrom}.sorted.vcf.gz.tbi
        } || { # except
            echo "tabix failed, initiate sorting/tabixing..."
            bcftools sort ${fp} -T ${TEMP} -Oz -o ./${inprefix}.chr${chrom}.sorted.vcf.gz
            tabix -p vcf ./${inprefix}.chr${chrom}.sorted.vcf.gz
        }
    else
        echo "input tabixed/sorted, skip tabixing/sorting."
        cp ${fp} ./${inprefix}.chr${chrom}.sorted.vcf.gz
        cp ${fp}.tbi ./${inprefix}.chr${chrom}.sorted.vcf.gz.tbi
    fi
}
export -f CHECK_SORT


SPLIT_CHR () {
    bcftools view ${1} -r ${2} -Oz -o ${1/.vcf.gz/}.chr${2}.sorted.vcf.gz
}
export -f SPLIT_CHR


LIFT_OVER () {
    lift=$PBS_O_WORKDIR/required_tools/lift/LiftMap.py
    cpath=$PBS_O_WORKDIR/required_tools/chainfiles
    cfilename=$(grep "Use chain file" $buildcheck | tr -d ' ' | tr ':' '\t' | tr -d '"' | cut -f 2 | sed -e 's/->/ /g')
    nchains=$(echo $cfilename | tr ' ' '\n' | wc -l | awk '{print $1}')
    checknone=$(grep "Use chain file" $buildcheck | grep "none" | wc -l)

    name=${inprefix}.chr$1
        
    echo "Remove multi-allelic variants"
    bcftools view $name.sorted.vcf.gz -M 2 -m 2 | bcftools norm /dev/stdin -d both -Oz -o $name.sorted.bi.vcf.gz
    tabix -p vcf $name.sorted.bi.vcf.gz

    echo "Converting vcf to ped/map.."
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' $name.sorted.bi.vcf.gz > $name.sorted.bi.pos
    # plink --vcf $name.sorted.bi.vcf --make-bed --a1-allele $name.sorted.bi.pos 5 3 --biallelic-only strict --set-missing-var-ids @:#:\$1:\$2 --vcf-half-call missing --out $name.sorted.bi
    # plink --bfile $name.sorted.bi --recode --a1-allele $name.sorted.bi.bim 5 2 --double-id --set-missing-var-ids @:#:\$1:\$2 --out $name
    plink --vcf $name.sorted.bi.vcf.gz --make-bed --a1-allele $name.sorted.bi.pos 5 3 --biallelic-only strict --set-missing-var-ids @:#:\$1:\$2 --vcf-half-call missing --double-id --recode ped --id-delim '_' --out $name
    # NOTE: '_' is FID_IID deliminator, keep watching if exception ID appeared.
    # TODO update to plink2 when it supported ped files.
    
    if [ $checknone -eq 1 ]; then
        echo "The data set is already based on the correct reference build (Grch37). Just converting and copying it."
        # plink --bfile $name.sorted.bi --make-bed --a1-allele $name.sorted.bi.bim 5 2 --double-id --set-missing-var-ids @:#:\$1:\$2 --allow-extra-chr --out $name.lifted_already_GRCh37.sorted.with_dup
        plink2 --vcf $name.sorted.bi.vcf.gz --make-bed --double-id --set-missing-var-ids @:#:\$1:\$2 --allow-extra-chr --out $name.lifted_already_GRCh37.sorted.with_dup
        
        cut -f 2 $name.lifted_already_GRCh37.sorted.with_dup.bim | sort | uniq -d > $name.list_multi_a_markers.txt
        ndup=$(wc -l $name.list_multi_a_markers.txt | awk '{print $1}')
        
        if [ "$ndup" -ge 1 ]; then
            echo "Still found duplicate variant ids or multiallelic markers after filtering. Performing additional filtering."
            dupflag=$(echo -e "--exclude $name.list_multi_a_markers.txt")
        else
            dupflag=""
        fi
        # plink --bfile $name.lifted_already_GRCh37.sorted.with_dup $dupflag --a1-allele $name.sorted.bi.bim 5 2 --make-bed --out $name.lifted_already_GRCh37
        plink2 --bfile $name.lifted_already_GRCh37.sorted.with_dup $dupflag --make-bed --out $name.lifted_already_GRCh37
        
        for bfile in bed bim fam; do
            mv $name.lifted_already_GRCh37.${bfile} $inprefix.lifted_already_GRCh37.chr$1.${bfile}
        done

        return 1
    fi

    echo "Number of chain files: $nchains"
    echo "Lifting ped/map using chain file(s): $cfilename..."
    if [ "$nchains" -eq 1 ]; then
         liftname=$(echo $cfilename |  sed -e 's/\.chain.*//g')
         $lift -p $name.ped -m $name.map -c $cpath/$cfilename -o $name.lifted_$liftname

    fi
    if [ "$nchains" -eq 2 ]; then
         first=$(echo $cfilename | tr ' ' '\n' | head -n 1)
         second=$(echo $cfilename | tr ' ' '\n' | tail -n 1)
         liftname=$(echo $first |  sed -e 's/\.chain.*//g')
         $lift -m $name.map -c $cpath/$first -o $name.lifted_$liftname
         liftname=$(echo $second |  sed -e 's/\.chain.*//g')
         $lift -m $name.lifted_$liftname.map -c $cpath/$second -o $name.lifted_$liftname
    fi

    echo "Converting lifted output $name.lifted_$liftname to bed/bim/fam.."
    plink --file $name.lifted_$liftname --a1-allele $name.sorted.bi.pos 5 3 --double-id --set-missing-var-ids @:#:\$1:\$2 --allow-extra-chr --make-bed --out $name.lifted_$liftname.sorted.with_dup
    # TODO: update to plink2 when it supported ped files.   

    cut -f 2 $name.lifted_$liftname.sorted.with_dup.bim | sort | uniq -d > $name.list_multi_a_markers.txt
        ndup=$(wc -l $name.list_multi_a_markers.txt | awk '{print $1}')
        if [ "$ndup" -ge 1 ]; then
            echo "Still found duplicate variant ids or multiallelic markers after filtering. Performing additional filtering."
            dupflag=$(echo -e "--exclude $name.list_multi_a_markers.txt")
            # plink --bfile $name.lifted_$liftname.sorted.with_dup $dupflag --make-bed --a1-allele $name.sorted.bi.bim 5 2 --allow-extra-chr --out $name.lifted_$liftname
            plink2 --bfile $name.lifted_$liftname.sorted.with_dup $dupflag --make-bed --allow-extra-chr --out $name.lifted_$liftname
        else
            echo "No duplicate. Skip additional plink, move the filename directly."
            dupflag=""
            # plink --bfile $name.lifted_$liftname.sorted.with_dup $dupflag --make-bed --a1-allele $name.sorted.bi.bim 5 2 --allow-extra-chr --out $name.lifted_$liftname
            plink2 --bfile $name.lifted_$liftname.sorted.with_dup $dupflag --make-bed --allow-extra-chr --out $name.lifted_$liftname
            for bfile in bed bim fam; do
                mv $name.lifted_$liftname.sorted.with_dup.${bfile} $inprefix.lifted_$liftname.chr$1.${bfile}
            done
        fi
    
    rm $name.map
    rm *with_dup*
    rm $name.lifted_$liftname.sorted.*
    rm $name.list_multi_a_markers.txt
    rm $name.sorted.bi.*
}
export -f LIFT_OVER


if [ ! -z $custom_temp ]; then
    mkdir -p $custom_temp/$PBS_JOBID
    export TEMP=$custom_temp/$PBS_JOBID
else
    mkdir -p $PBSTMPDIR/$PBS_JOBID
    export TEMP=$PBSTMPDIR/$PBS_JOBID
fi

if [ ! -d $myoutdir ]; then
    mkdir -p $myoutdir
fi

cd $myoutdir


vcfstarttime=$(date +%s)
# Test if input came as split chromosome
if [[ ${myinput} == *.txt ]]; then
    echo "Input autosome list detected, run parallel pipeline."
    parallel CHECK_SORT {1} :::: ${myinput}
elif [[ ${myinput} == *.vcf.gz ]]; then 
    echo "Input file is merged; split by chromsome."
    CHECK_SORT "ALL\t${myinput}"
    SPLIT_CHR ./${inprefix}.chrALL.sorted.vcf.gz {1} ::: {1..22}
else  
    echo "Invalid file format, please check the input."
    exit
fi
vcfendtime=$(date +%s)


liftstarttime=$(date +%s)
parallel LIFT_OVER {1} ::: {1..22}
liftendtime=$(date +%s)


vcfruntime=$((vcfendtime-vcfstarttime))
liftruntime=$((liftendtime-liftstarttime))
echo "vcf convert time: $vcfruntime"
echo "liftover run time: $liftruntime"

#output example
#$inprefix.lifted_$liftname.chr${chrom}.${bfile}