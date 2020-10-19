module load plink
module load R/3.6.1

cd /scratch/midway2/liliw1/

# dataset DGN: dir_plink=/project2/xuanyao/data/DGN
# dataset DGN: dir_txt=/project2/xuanyao/data/DGN/txt_file_of_DGN
# dataset DGN: dir_script=/project2/xuanyao/data/DGN
# dataset DGN: awk 'BEGIN{print "id\tchr\tpos"}{printf("%s\tchr%s\t%s\n",$2,$1,$4)}'

dir_plink=/project2/xuanyao/xuanyao/cancer_gbat/breast_cancer
dir_txt=/project2/xuanyao/llw/breastcancerTCGA/txt_file
dir_script=/scratch/midway2/liliw1/TCGA


for CHR in `seq 22`
do
# convert bfiles to .raw files
plink \
--bfile $dir_plink/genotype/chr${CHR}_QCed \
--recodeA \
--out $dir_txt/chr${CHR}

# extract genotype meta
awk \
'BEGIN{print "id\tchr\tpos"}{key=sprintf("%s:%s", $1, $4);printf("%s\tchr%s\t%s\n",key,$1,$4)}' \
$dir_plink/genotype/chr${CHR}_QCed.bim > \
$dir_txt/chr${CHR}_genotype_meta.txt

# convert genotype to txt file
Rscript --no-restore --no-save $dir_script/make.geno.R $CHR

done
