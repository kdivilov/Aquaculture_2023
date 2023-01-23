library(ggbiplot)

genes_dCq = read.csv("genes_dCq.csv")

#Biplot
pca = prcomp(t(scale(t(genes_dCq[,4:10]))), scale = T)
ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = factor(genes_dCq$Chr8.9719736)) +
  scale_color_discrete(name = 'Chr8:9719736', labels=c("ref/ref","ref/alt","alt/alt"))

#Subset of families with >2 homozygous (ref) and >2 heterozygous individuals
genes_dCq_subset = genes_dCq[which(genes_dCq$Family %in% c(29.001,29.016,29.017,29.025,29.042,29.048,29.056,29.072,29.073,29.077)),]

#Fold change in gene expression between individuals heterozygous at the Chr8:9719736 SNP and those homozygous (ref) at the Chr8:9719736 SNP 
genes_FC = data.frame(matrix(NA,10,7))
rownames(genes_FC) = unique(genes_dCq_subset$Family)
names(genes_FC) = c("ADAR","IFI44","IRF2","IRF8","MDA5","STING","Viperin")

for(i in 1:nrow(genes_FC)){
  for(j in 1:ncol(genes_FC)){
    mean_heterozygous = mean(genes_dCq_subset[which(genes_dCq_subset$Family==rownames(genes_FC)[i] & genes_dCq_subset$Chr8.9719736==1),names(genes_FC)[j]])
    mean_homozygous = mean(genes_dCq_subset[which(genes_dCq_subset$Family==rownames(genes_FC)[i] & genes_dCq_subset$Chr8.9719736==0),names(genes_FC)[j]])
    genes_FC[i,j] = mean_heterozygous/mean_homozygous
  }
}

#Mean fold change in gene expression
colMeans(genes_FC)