module load plink

for CHR in `seq 22`
do
# convert bfiles to .raw files
plink \
--bfile /project2/xuanyao/data/DGN/genotype/chr${CHR}_QCed \
--recodeA \
--out /project2/xuanyao/data/DGN/txt_file_of_DGN/chr${CHR}

# extract genotype meta
awk \
'BEGIN{print "id\tchr\tpos"}{printf("%s\tchr%s\t%s\n",$2,$1,$4)}' \
/project2/xuanyao/data/DGN/genotype/chr${CHR}_QCed.bim > \
/project2/xuanyao/data/DGN/txt_file_of_DGN/chr${CHR}_genotype_meta.txt

# convert genotype to txt file
Rscript --no-restore --no-save /project2/xuanyao/data/DGN/make_geno.R $CHR
done



