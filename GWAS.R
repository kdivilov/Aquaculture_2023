library(vcfR)

vcf = read.vcfR("mbp_snps_imputed.vcf.gz")
geno = extract.gt(vcf, element = "GT")
colnames(geno) = gsub("\\.bam*","",basename(colnames(geno)))
geno[geno=="0/0"] = 0
geno[geno=="0/1"] = 1
geno[geno=="1/1"] = 2
geno = t(geno)
mode(geno) = "numeric"

pedigree = data.frame(names(which(table(gsub("_.*","",rownames(geno)))==2)))
names(pedigree) = "Family"
pedigree$Female = paste(pedigree$Family,"F",sep="_")
pedigree$Male = paste(pedigree$Family,"M",sep="_")

simulate_cross = function(x){
  if(x[1]==0 && x[2]==0){
    return(0)
  }
  if(x[1]==2 && x[2]==2){
    return(2)
  }
  if(x[1]==0 && x[2]==2){
    return(1)
  }
  if(x[1]==2 && x[2]==0){
    return(1)
  }
  if(x[1]==1 && x[2]==1){
    return(1)
  }
  if(x[1]==0 && x[2]==1){
    return(0.5)
  }
  if(x[1]==1 && x[2]==0){
    return(0.5)
  }
  if(x[1]==2 && x[2]==1){
    return(1.5)
  }
  if(x[1]==1 && x[2]==2){
    return(1.5)
  }
}

geno_sim = matrix(NA,nrow(pedigree),ncol(geno))
rownames(geno_sim) = pedigree$Family
colnames(geno_sim) = colnames(geno)
for(i in 1:nrow(geno_sim)){
  geno_sim[i,] = apply(geno[match(c(pedigree$Female[i],pedigree$Male[i]),rownames(geno)),],2,simulate_cross)
}

pheno = read.csv("phenotypes.csv")
pheno$Family = factor(pheno$Family)
pheno$Year = factor(pheno$Year)

pheno_mean = aggregate(Status~Family,pheno,mean)
pheno_mean$Year = NA
pheno_mean$Year[grep("27.",pheno_mean$Family)] = "2018"
pheno_mean$Year[grep("29.",pheno_mean$Family)] = "2020"
pheno_mean = pheno_mean[-which(is.na(pheno_mean$Year)),]
pheno_mean = pheno_mean[which(pheno_mean$Family %in% pedigree$Family),]

gemma_pheno = pheno_mean$Status*100
gemma_covar = model.matrix(~as.factor(pheno_mean$Year))
gemma_geno = t(geno_sim)
gemma_geno = cbind(rownames(gemma_geno),vcf@fix[,5],vcf@fix[,4],gemma_geno)
gemma_anno = cbind(rownames(gemma_geno),
                   sapply(strsplit(rownames(gemma_geno),"_"), `[`, 2),
                   sapply(strsplit(sapply(strsplit(rownames(gemma_geno),"_"), `[`, 1),"Chr"), `[`, 2))

write.table(gemma_pheno,"gemma_pheno.txt",row.names = F,col.names = F)
write.table(gemma_covar,"gemma_covar.txt",row.names = F,col.names = F)
write.table(gemma_geno,"gemma_geno.txt",row.names = F,col.names = F,quote = F)
write.table(gemma_anno,"gemma_anno.txt",row.names = F,col.names = F,quote = F)

#GWAS (run in command line)
./gemma-0.98.5 -g gemma_geno.txt -p gemma_pheno.txt -c gemma_covar.txt -a gemma_anno.txt -gk -maf 0.05 -o gemma_K
./gemma-0.98.5 -g gemma_geno.txt -p gemma_pheno.txt -c gemma_covar.txt -a gemma_anno.txt -k output/gemma_K.cXX.txt -maf 0.05 -lmm 2 -o gemma_gwas