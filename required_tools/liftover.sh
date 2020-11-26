#!/bin/bash

# Example: bash liftover.sh test.vcf.gz

myinput=$1
inprefix=${myinput/.vcf.gz/}

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${myinput} > ${inprefix}.pos
plink --vcf ${myinput} --make-bed --a1-allele ${inprefix}.pos 5 3 --biallelic-only strict --set-missing-var-ids @:#:\$1:\$2 --vcf-half-call missing --double-id --recode ped --out ${inprefix}

lift/LiftMap.py -p ${inprefix}.ped -m ${inprefix}.map -c chainfiles/GRCh37_to_GRCh38.chain -o ${inprefix}.lifted_GRCh37_to_GRCh38

plink --file ${inprefix}.lifted_GRCh37_to_GRCh38 --a1-allele ${inprefix}.pos 5 3 --double-id --set-missing-var-ids @:#:\$1:\$2 --allow-extra-chr --make-bed --out ${inprefix}.lifted_GRCh37_to_GRCh38.sorted.with_dup


