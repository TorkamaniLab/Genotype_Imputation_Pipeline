#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=540:00:00
#SBATCH --mem=100G


#if running outside pipeline (standalone version), create symbolic links to make your file names match to the expected syntax
#the file suffix should look like this: .lifted_hg38_to_hg19.GH.chr1.vcf.gz
#chromosome number should be the only variable part of the name
#Example from my work dir:
#wdir=/mnt/stsi/stsi5/raqueld/HGDP/WGS/hg19_from_hg38
#cd $wdir/ancestry
#for i in $wdir/picard_QC/hgdp_wgs.20190516.statphase.autosomes.chr{1..22}.hg19.concat.rename.bi.dense.vcf.gz; do
#    chr=$(basename $i | sed -e 's/.*\.chr//g' | sed -e 's/\..*//g')
#    newname=$(basename $i | sed -e "s/chr$chr\.//g" | sed -e "s/\.vcf\.gz$/\.lifted_hg38_to_hg19\.GH\.chr$chr\.vcf\.gz/g")
#    ln -s $i $newname
#    ln -s $i.tbi $newname.tbi
#done
#then run the job
#sbatch --export=myinput=$wdir/ancestry/hgdp_wgs.20190516.statphase.autosomes.hg19.concat.rename.bi.dense.lifted_hg38_to_hg19.GH,myoutdir=$wdir/ancestry,WGS=yes --job-name=ancestry_HGDP 3_ancestry_analysis.slurm.sh

#Step 3 consists of
#Prune markers in input file by 0.05 LD threshold
#Intersect input pruned markers with reference panel (1KG)
#Perform ancestry analysis
#Split subjects by ancestry (threshold >=0.95)


date
echo "Running on node:"
hostname
pwd

echo SLURM_SUBMIT_DIR: $SLURM_SUBMIT_DIR

module purge
module load vcftools
module load admixture
module load samtools
module load R


export plink="$SLURM_SUBMIT_DIR/required_tools/plink"
export plink2="$SLURM_SUBMIT_DIR/required_tools/plink2"
export split_by_ancestry="$SLURM_SUBMIT_DIR/required_tools/split_by_ancestry/split_by_ancestry.R"


starttime=$(date +%s)

export filename=$(basename $myinput)
export inprefix=$(basename $myinput | sed -e 's/\.vcf.gz$//g')
export outsubdir=$(basename $myinput | sed -e 's~\.lifted.*~~g')
export indir=$(dirname $myinput)

if [ ! -d $myoutdir/$outsubdir ]; then
    mkdir -p $myoutdir/$outsubdir
fi

cd $myoutdir/$outsubdir

echo "Preparing data to run ancestry analysis."

vcfstarttime=$(date +%s)


PRUNE_FUN() {
    echo "Pruning markers using 0.05 LD threshold"
    if [ ${WGS} == "yes" ]; then
        echo "Extracting 23andMe positions only for ancestry analysis, for speeding up ancestry analysis. Original positions will be restored later"
        #positions extracted from 190410_snps.23andme.clean.drop_dup.sorted.data
        bcftools view -R $SLURM_SUBMIT_DIR/required_tools/chr_pos_23andme.txt ${myinput}.chr$1.vcf.gz -Ov -o ${inprefix}.chr$1.23andMe_pos.vcf
        $plink2 --vcf ${inprefix}.chr$1.23andMe_pos.vcf --indep-pairwise 100 10 0.05 --out $inprefix.chr$1 # produces <name.prune.in> and <name.prune.out>
    else
        $plink2 --vcf $myinput.chr$1.vcf.gz --indep-pairwise 100 10 0.05 --out $inprefix.chr$1 # produces <name.prune.in> and <name.prune.out>
    fi

    # plink --bfile $indir/$inprefix --extract $inprefix.prune.in --recode vcf bgz --keep-allele-order --set-missing-var-ids @:#\$1:\$2 --out $inprefix.pruned
    vcftools --gzvcf $myinput.chr$1.vcf.gz --snps $inprefix.chr$1.prune.in --recode --recode-INFO-all --out $inprefix.chr$1.pruned

    bgzip -c $inprefix.chr$1.pruned.recode.vcf > $inprefix.chr$1.pruned.vcf.gz
    tabix -p vcf $inprefix.chr$1.pruned.vcf.gz
    
    rm $inprefix.chr$1.pruned.recode.vcf
}
export -f PRUNE_FUN

echo "Running pruning..."
parallel PRUNE_FUN {1} ::: {1..22}
echo "pruning done."


ISEC_FUN() {
    echo "Generating intersection between input data and reference panel for ancestry"
    bcftools isec \
      $inprefix.chr$1.pruned.vcf.gz \
      /mnt/stsi/stsi0/raqueld/1000G/ALL.chr$1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.clean.vcf.gz \
      -p ${inprefix}.chr$1_tmp -n =2  -w 1,2 -Oz
    
    bcftools merge ${inprefix}.chr$1_tmp/0000.vcf.gz ${inprefix}.chr$1_tmp/0001.vcf.gz \
      -Oz -o $inprefix.chr$1.pruned.intersect1KG.vcf.gz
    tabix -f -p vcf $inprefix.chr$1.pruned.intersect1KG.vcf.gz
}
export -f ISEC_FUN

echo "Running isec..."
parallel ISEC_FUN {1} ::: {1..22}
echo "isec done."


echo "Concatenate intersection for glocal ancestry inference"
ls -v $inprefix.chr*.pruned.intersect1KG.vcf.gz > $inprefix.pruned.intersect1KG.vcflist
bcftools concat -f $inprefix.pruned.intersect1KG.vcflist -Ov -o $inprefix.pruned.intersect1KG.vcf

# rename 1000G sample ID as FID_IID, required by the make-bed later
echo "Rename the reference panel sample ID"
bcftools query -l ${inprefix}.chr1_tmp/0001.vcf.gz > 0001.txt
while read line; do printf "${line}\t${line}_${line}\n"; done < 0001.txt > 0001.rename
bcftools reheader $inprefix.pruned.intersect1KG.vcf -s 0001.rename > $inprefix.pruned.intersect1KG.id-delim.vcf
rm $inprefix.pruned.intersect1KG.vcf

echo "Converting intersection to bed format"
# plink --vcf $inprefix.pruned.intersect1KG.vcf.gz --keep-allele-order --make-bed --out $inprefix.pruned.intersect1KG
$plink2 --vcf $inprefix.pruned.intersect1KG.id-delim.vcf --make-bed --id-delim '_' --out $inprefix.pruned.intersect1KG
rm inprefix.pruned.intersect1KG.id-delim.vcf


echo "Generating population file for admixture"
inputN=$(bcftools query -l $inprefix.chr22.pruned.vcf.gz | wc -l)

if [ -f $inprefix.pruned.intersect1KG.pop ]; then
      rm $inprefix.pruned.intersect1KG.pop
      echo
fi

for i in $(seq 1 1 $inputN); do
     echo
done >> $inprefix.pruned.intersect1KG.pop

cat /gpfs/group/torkamani/raqueld/ancestry/1000G_P3_super_pop.pop  >> ./$inprefix.pruned.intersect1KG.pop


vcfendtime=$(date +%s)


echo "Glocal ancestry inference by ADMIXTURE"
admstarttime=$(date +%s)
    admixture --supervised $inprefix.pruned.intersect1KG.bed 5 -j16
admendtime=$(date +%s)


poststarttime=$(date +%s)

echo "Get subject IDs with ancestry inference results, including both input and reference subjects"
cat $inprefix.pruned.intersect1KG.fam | awk '{print $1 "_" $2}' > $inprefix.pruned.intersect1KG.subjectIDs
paste $inprefix.pruned.intersect1KG.subjectIDs $inprefix.pruned.intersect1KG.5.Q | tr '\t' ' ' > $inprefix.pruned.intersect1KG.5.Q.IDs

echo "Get subject IDs"
# Shaun: retire the old path, input was split vcf.gz file now
# cat $inprefix.fam | awk '{print $1 "_" $2 "\t" $1 "\t" $2}' > $inprefix.pruned.subjectIDs
bcftools query -l $inprefix.chr1.pruned.vcf.gz | tr '_' '\t' | awk '{print$1"_"$2"\t"$1"\t"$2}' > $inprefix.pruned.subjectIDs

echo "Split by ansestry cutoff 0.95"
Rscript $split_by_ancestry $inprefix.pruned.intersect1KG.5.Q.IDs $inprefix.pruned.subjectIDs 0.95


SPLIT_FUN() {
    # bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' $myinput > $inprefix.pos
    # plink --vcf $myinput --a1-allele $inprefix.pos 5 3 --make-bed --out $inprefix
    $plink2 --vcf $myinput.chr$1.vcf.gz --make-bed --id-delim '_' --out $inprefix.chr$1

    for i in {1..5} mixed; do 
        if [ -f $inprefix.pruned.intersect1KG.5.Q.IDs.$i.ids ]; then
            # plink --bfile $inprefix --keep $inprefix.pruned.intersect1KG.5.Q.IDs.$i.ids --a1-allele $inprefix.pos 5 3 --make-bed --out $inprefix.ancestry-$i
            $plink2 --bfile $inprefix.chr$1 --keep $inprefix.pruned.intersect1KG.5.Q.IDs.$i.ids --make-bed --out $inprefix.ancestry-$i.chr$1
        fi
        echo "ancestry-${i}.chr$1 done"
    done

    # Shaun: mixed included in the for-loop above, allowed by bash-syntax
    # if [ -f $inprefix.pruned.intersect1KG.5.Q.IDs.mixed.ids ]; then
    #     # plink --bfile $inprefix --keep $inprefix.pruned.intersect1KG.5.Q.IDs.mixed.ids --a1-allele $inprefix.pos 5 3 --make-bed --out $inprefix.ancestry-mixed
    #     $plink2 --bfile $inprefix.chr$1 --keep $inprefix.pruned.intersect1KG.5.Q.IDs.mixed.ids --make-bed --out $inprefix.ancestry-mixed.chr$1
    #     echo "ancestry-mixed.chr$1 done"
    # fi
}
export -f SPLIT_FUN

echo "Converting input to bed format"
parallel SPLIT_FUN {1} ::: {1..22}


postendtime=$(date +%s)


echo "Execution complete. Results ready for HWE filtering and phasing at $myoutdir/$outsubdir/$inprefix.ancestry-_.bed"


vcfruntime=$((vcfendtime-vcfstarttime))
admruntime=$((admendtime-admstarttime))
postruntime=$((postendtime-poststarttime))

echo "Preprocessing run time: $vcfruntime"
echo "Admixture run time: $admruntime"
echo "Postprocessing run time: $postruntime"
