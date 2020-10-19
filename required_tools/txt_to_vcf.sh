
echo '##fileformat=VCFv4.2'
echo '##FILTER=<ID=PASS,Description="All filters passed">'
echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
echo "##bash make_vcf.sh ${1}; Date=$(date '+%d/%m/%Y %H:%M:%S')"
printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMASKED_SAMPLE\n"

while read line; do
	chrom=$(echo $1 | awk '{print$1}')
	bp=$(echo $1 | awk '{print$2}')
	printf "${chrom}\t${bp}\t.\t"
chr22   15382411        .       G       A       255     SVM     .       GT      0/0
done < $1
