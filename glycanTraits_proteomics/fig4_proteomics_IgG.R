## Tiphaine Martin - 6/6/2019
## Script developped for the paper doi: https://doi.org/10.1101/662957

## if you use this script as template, please cite this paper
## doi: https://doi.org/10.1101/662957

############################# 
#        Figure 4
#############################
##Correlation between proteins associated with IgG
igG_proteins <-read.table(file="/figures_papers/data/IgG_proteins.csv",header=TRUE,sep=",")
dim(igG_proteins)
rownames(igG_proteins) <- igG_proteins$X
igG_proteins<- igG_proteins[,-1]
igG_proteins <- as.matrix(igG_proteins)
coul=col3(100)
#coul = colorRampPalette(brewer.pal(8, "RdBu"))(25)
heatmap(igG_proteins,col= coul)
colRwB=colorRamp2(c(-0.05, -0.04,-0.03,-0.02,-0.01,
                    0, 
                    0.01,0.02,0.03,0.04,0.05), 
                  c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                    "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                    "#D6604D", "#B2182B", "#67001F"))
Heatmap(igG_proteins,col=colRwB,color_space = colRwB)

## Correlation betweeen proteins
#/anaconda3/lib/python3.6/site-packages
use_python("/anaconda3/bin/python3.6")
devtools::install_github("ropenscilabs/umapr")
library(umapr)
library(tidyverse)

# select only numeric columns for proteins data
df <- AITD_prot_sub3[,20:1148]
pdf <- t(df)
# run UMAP algorithm
pembedding <- umap(pdf)
# plot the result
ggplot(pembedding, aes(UMAP1, UMAP2)) + geom_point() +theme_bw()

## correlation of TSH between two assays
## all
fit <- lm(SL000589 ~ TSH, data = AITD_prot_sub)
summary(fit)

table(AITD_prot_sub$batch)
AITD_prot_subR <- AITD_prot_sub[which(AITD_prot_sub$batch %in% "Roche"),]
fit <- lm(SL000589 ~ TSH, data = AITD_prot_subR)
summary(fit)

AITD_prot_subA <- AITD_prot_sub[which(AITD_prot_sub$batch %in% "Abbott"),]
fit <- lm(SL000589 ~ TSH, data = AITD_prot_subA)
summary(fit)

fit1 <- lm(log10(SL000589) ~ log10(TSH), data = AITD_prot_sub)
summary(fit1)

fit1 <- lm(log10(SL000589) ~ log10(TSH), data = AITD_prot_subR)
summary(fit1)
fit1 <- lm(log10(SL000589) ~ log10(TSH), data = AITD_prot_subA)
summary(fit1)

p11 <- ggplot(AITD_prot_sub, aes(x=SL000589, y=TSH,colour = batch, shape = batch)) + 
  geom_point() + geom_smooth(method=lm) +
  scale_x_log10() + scale_y_log10() +
  annotate("text", x=3000, y=0.75, label = "Roche: R^2=0.43") + 
  annotate("text", x=3000, y=0.5, label = "Abbott: R^2=0.47") + 
  theme_bw()
p11

