# Genotype_Imputation_Pipeline
> Final version of the pipeline for processing TOPMED and other dbGaP datasets, includes: 
automatic identification of reference genome build and selection of correct chain file, mapping to GRCh37 reference build, initial quality control, ancestry analysis, post-ancestry quality control, phasing, imputation, post-imputation quality control.

## Dependencies

__This pipeline ONLY on Garibaldi__  
```module load vcftools```  
```module load plink2/1.90b3.42```  
```module load samtools/1.9```  
```module load R```  
```~/split_by_ancestry/split_by_ancestry.R``` (copied /gpbs/home/raqueld/split_by_ancestry)  
```~/bin/LiftMap.py``` (copied /gpfs/home/raqueld/bin)  

## Usage example

### Step 0
automatic identification of reference genome build and selection of correct chain file  
```ruby
qsub 0_check_vcf_build.job -v  myinput=/stsi/raqueld/vcf/6800_JHS_all_chr_sampleID_c2.vcf,myoutput=/stsi/raqueld/0_check_vcf_build/6800_JHS_all_chr_sampleID_c2.BuildChecked,copyoutput=yes,gz=yes -N 0_6800_JHS_all_chr_sampleID_c2```

- **Option 1**
    - 
    
### Step 1
mapping to GRCh37 reference build  
```qsub 1_lift_vcfs_to_GRCh37.job -v myinput=/stsi/raqueld/vcf/6800_JHS_all_chr_sampleID_c2.vcf,buildcheck=/stsi/raqueld/0_check_vcf_build/6800_JHS_all_chr_sampleID_c2.BuildChecked,myoutdir=/stsi/raqueld/1_lift,copyoutput=yes -N 1_6800_JHS_all_chr_sampleID_c2```


### Step 2
initial quality control  
```qsub 2_Genotype_Harmonizer.job -v myinput=/gpfs/home/raqueld/mapping_MESA/mesa_genotypes-black.lifted_NCBI36_to_GRCh37.bed,myoutdir=/gpfs/home/raqueld/mapping_MESA -N 2_N_GH.mesa_genotypes-black```

### Step 3
ancestry analysis  
```qsub 3_ancestry_analysis.job -v myinput=/stsi/raqueld/2_GH/6800_JHS_all_chr_sampleID_c1.lifted_hg19_to_GRCh37.GH.fix.vcf.gz,myoutdir=/stsi/raqueld/3_ancestry -N 3_6800_JHS_all_chr_sampleID_c1```

### Step 4
post-ancestry quality control  
```qsub 4_split_QC2.job -v myinput=/gpfs/home/raqueld/mapping_MESA/mesa_genotypes-black.lifted_NCBI36_to_GRCh37.GH.bed,myoutdir=/stsi/raqueld/N_tests,hwe='',geno=0.1,mind=0.1 -N 4_N_mesa_genotypes-black```
```qsub 4_split_QC2.job -v myinput=/stsi/raqueld/3_ancestry/6800_JHS_all_chr_sampleID_c1/6800_JHS_all_chr_sampleID_c1.lifted_hg19_to_GRCh37.GH.ancestry-5.bed,myoutdir=/stsi/raqueld/4_split_QC2,hwe='',geno=0.1,mind=0.1 -N 4_6800_JHS_all_chr_sampleID_c1```

### Step 5
phasing  
```qsub 5_phase.job -v myinput=/stsi/raqueld/N_tests/aric_genotypes-black/aric_genotypes-black.lifted_NCBI36_to_GRCh37.GH.chr1.bed,myoutdir=/stsi/raqueld/5_N_tests,reftype=HRC -N 5_N_mesa_genotypes-black```

  - the imput must have the suffix *.lifted*.chr1.bed, *lifted*.chr2.bed, *.lifted*.chr3.bed, etc. The previous steps in the pipeline generate those suffixes automatically, but keep these suffixes in mind if you are running this step as a stand alone tools, without running the previous steps

### Step 6
imputation, post-imputation quality control  
```qsub 6_impute.job -v myinput=/stsi/raqueld/5_N_tests/mesa_genotypes-white/mesa_genotypes-white.lifted_NCBI36_to_GRCh37.GH.chr18.phased.vcf.gz,myoutdir=/stsi/raqueld/6_N_tests,reftype=HRC -N 6_mesa_genotypes-white.lifted_NCBI36_to_GRCh37.GH.chr18```
  - if running this step as stand alone tool, input file must have the suffix .lifted_*.chr*.phased.vcf.gz, otherwise the pipeline wont work, if you use the previous step to generate this input file, then it will work fine.
date


## Meta

Raquel Dias – [@RaquelDiasSRTI](https://twitter.com/RaquelDiasSRTI) – raqueld@scripps.edu
Shaun Chen - [@ShaunFChen](http://twitter.com/ShaunFChen) - sfchen@scripps.edu


## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request
