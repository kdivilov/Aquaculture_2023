mkdir mbp
mkdir mbp/raw
mkdir mbp/trim
mkdir mbp/map
mkdir mbp/vcf

#Cgigas_RefSeq.fasta is GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fasta (RefSeq GCF_902806645.1) from NCBI with linkage groups renamed to chromosomes according to Peñaloza et al. (2021) and unplaced scaffolds removed
bwa index Cgigas_RefSeq.fasta

#Download and place raw reads from NCBI BioProject PRJNA873124 into the mbp/raw folder

#Trim raw reads
sh mbp_trim.sh

#Map trimmed reads to the reference genome
sh mbp_map.sh

#Call biallelic SNPs
bcftools mpileup -f Cgigas_RefSeq.fasta mbp/map/*.bam -C 50 -Q 30 -q 40 -a "AD,DP" -I -Ou | bcftools call -vm -Ou | bcftools view -m2 -M2 -Oz -o mbp/vcf/mbp_snps.vcf.gz

#Impute missing SNP genotypes
java -jar LinkImputeR.jar -s LinkImputeR_config.ini
java -jar LinkImputeR.jar mbp/vcf/LinkImputeR.xml "Case 1" mbp/vcf/mbp_snps_imputed.vcf.gz

